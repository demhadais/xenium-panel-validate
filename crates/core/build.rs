use std::{
    collections::HashSet,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    str::FromStr,
};

use anyhow::Context;
use bstr::ByteSlice;
use noodles::gff::feature::{RecordBuf, record_buf::attributes::field::Value};
use serde::Deserialize;
use url::Url;

// About 60,000 genes for human (the larger set of annotations), the closest power of 2 is 2^16
const N_GENES: usize = 65_356;

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    if !option_env!("BUILD_XENIUM_PANEL_VALIDATE")
        .map(bool::from_str)
        .transpose()?
        .is_some_and(|build| build)
    {
        return Ok(());
    }

    let Config { human, mouse } = toml::from_slice(include_bytes!("genes.toml"))
        .context("failed to parse config from genes.toml")?;

    let mut annotations = HashSet::with_capacity(N_GENES);
    read_gene_annotations_into(&human.gene_annotations_path, &mut annotations)?;

    let http_client = reqwest::Client::new();
    let unavailable_gene_sets = [
        human.xenium_v1_unavailable_genes_url,
        human.xenium_prime_unavailable_genes_url,
        mouse.xenium_v1_unavailable_genes_url,
        mouse.xenium_prime_unavailable_genes_url,
    ]
    .map(|url| fetch_unavailable_ensembl_ids(&http_client, url));
    let unavailable_gene_sets = futures::future::try_join_all(unavailable_gene_sets).await?;
    let [
        human_v1_unavailable,
        human_prime_unavailable,
        mouse_v1_unavailable,
        mouse_prime_unavailable,
    ] = unavailable_gene_sets.as_array().unwrap();

    let human_v1_map = construct_map(&annotations, human_v1_unavailable);
    write_map_to_file(
        &PathBuf::from("src/gene_list/chemistry/xenium_v1_human.rs"),
        "XENIUM_V1_HUMAN_GENES",
        &human_v1_map,
    )?;

    let human_prime_map = construct_map(&annotations, human_prime_unavailable);
    write_map_to_file(
        &PathBuf::from("src/gene_list/chemistry/xenium_prime_human.rs"),
        "XENIUM_PRIME_HUMAN_GENES",
        &human_prime_map,
    )?;

    annotations.clear();
    read_gene_annotations_into(&mouse.gene_annotations_path, &mut annotations)?;

    let mouse_v1_map = construct_map(&annotations, mouse_v1_unavailable);
    write_map_to_file(
        &PathBuf::from("src/gene_list/chemistry/xenium_v1_mouse.rs"),
        "XENIUM_V1_MOUSE_GENES",
        &mouse_v1_map,
    )?;

    let mouse_prime_enums = construct_map(&annotations, mouse_prime_unavailable);
    write_map_to_file(
        &PathBuf::from("src/gene_list/chemistry/xenium_prime_mouse.rs"),
        "XENIUM_PRIME_MOUSE_GENES",
        &mouse_prime_enums,
    )?;

    Ok(())
}

#[derive(Clone, Debug, Deserialize)]
struct Config {
    human: SpeciesConfig,
    mouse: SpeciesConfig,
}

#[derive(Clone, Debug, Deserialize)]
struct SpeciesConfig {
    xenium_v1_unavailable_genes_url: Url,
    xenium_prime_unavailable_genes_url: Url,
    gene_annotations_path: PathBuf,
}

type GeneAnnotationFileReader = noodles::gtf::io::Reader<BufReader<File>>;

fn read_gene_annotations_into(
    path: &Path,
    genes_buf: &mut HashSet<(String, String)>,
) -> anyhow::Result<()> {
    let mut reader =
        GeneAnnotationFileReader::new(BufReader::new(File::open(path).with_context(|| {
            format!(
                "failed to read gene annotations from {}",
                path.to_str().expect("path should be UTF-8")
            )
        })?));

    for record in reader.record_bufs() {
        let record = record?;

        let Some((ensembl_id, gene_name)) = parse_ensembl_id_and_name_from_gtf_record(&record)
        else {
            continue;
        };

        genes_buf.insert((ensembl_id, gene_name));
    }

    Ok(())
}

fn construct_map<'a>(
    genes: &'a HashSet<(String, String)>,
    unavailable_gene_ids: &HashSet<String>,
) -> phf_codegen::Map<'a, &'a str> {
    let mut map = phf_codegen::Map::new();

    for (ensembl_id, gene_name) in genes {
        if unavailable_gene_ids.contains(ensembl_id) {
            continue;
        }

        map.entry(ensembl_id.as_ref(), format!(r#""{gene_name}""#));
    }

    map
}

fn write_map_to_file(
    path: &Path,
    map_name: &str,
    map: &phf_codegen::Map<'_, &str>,
) -> anyhow::Result<()> {
    let file = File::create(path)
        .with_context(|| format!("failed to write file {}", path.to_str().unwrap()))?;
    let mut file_writer = BufWriter::new(file);

    writeln!(
        file_writer,
        "pub static {map_name}: phf::Map<&'static str, &'static str> = {};",
        map.build()
    )?;

    Ok(())
}

async fn fetch_unavailable_ensembl_ids(
    client: &reqwest::Client,
    url: Url,
) -> anyhow::Result<HashSet<String>> {
    #[derive(Deserialize)]
    struct Gene {
        gene_id: String,
    }

    let response = client.get(url).send().await?;

    let raw = response.bytes().await?;
    let mut reader = csv::Reader::from_reader(raw.iter().as_slice());

    let mut gene_ids = HashSet::with_capacity(1500);
    for row in reader.deserialize() {
        let Gene { gene_id } = row?;
        gene_ids.insert(gene_id);
    }

    Ok(gene_ids)
}

fn parse_ensembl_id_and_name_from_gtf_record(record: &RecordBuf) -> Option<(String, String)> {
    let attributes = record.attributes();

    let (Some(Value::String(gene_id)), Some(Value::String(gene_name))) =
        (attributes.get(b"gene_id"), attributes.get(b"gene_name"))
    else {
        unreachable!("'gene_id' and 'gene_name' should be strings");
    };

    Some((
        gene_id
            .to_str()
            .unwrap()
            .split('.')
            .next()
            .map(str::to_owned)
            .unwrap(),
        gene_name.to_str().map(str::to_owned).unwrap(),
    ))
}
