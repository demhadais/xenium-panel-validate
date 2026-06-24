use std::{
    collections::HashSet,
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
    str::FromStr,
};

use anyhow::Context;
use bstr::{BString, ByteSlice};
use heck::ToPascalCase;
use noodles::gff::feature::{RecordBuf, record_buf::attributes::field::Value};
use proc_macro2::Span;
use quote::quote;
use serde::Deserialize;
use syn::{Ident, LitStr};
use url::Url;

const N_GENES: usize = 24_000;

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    if !option_env!("BUILD_XENIUM_PANEL_VALIDATE")
        .map(bool::from_str)
        .transpose()?
        .is_some_and(|build| build)
    {
        return Ok(());
    }

    println!("cargo::rerun-if-changed=genes.toml");

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

    let human_v1_enums = make_enums("XeniumV1Human", &annotations, human_v1_unavailable);
    std::fs::write("src/gene_list/gene/xenium_v1_human.rs", human_v1_enums)?;

    let human_prime_enums = make_enums("XeniumPrimeHuman", &annotations, human_prime_unavailable);
    std::fs::write(
        "src/gene_list/gene/xenium_prime_human.rs",
        human_prime_enums,
    )?;

    annotations.clear();
    read_gene_annotations_into(&mouse.gene_annotations_path, &mut annotations)?;

    let mouse_v1_enums = make_enums("XeniumV1Mouse", &annotations, mouse_v1_unavailable);
    std::fs::write("src/gene_list/gene/xenium_v1_mouse.rs", mouse_v1_enums)?;

    let mouse_prime_enums = make_enums("XeniumPrimeMouse", &annotations, mouse_prime_unavailable);
    std::fs::write(
        "src/gene_list/gene/xenium_prime_mouse.rs",
        mouse_prime_enums,
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
    let mut reader = GeneAnnotationFileReader::new(BufReader::new(File::open(path)?));

    for record in reader.record_bufs() {
        let record = record?;

        let Some((ensembl_id, gene_name)) = parse_ensembl_id_and_name_from_gtf_record(&record)
        else {
            continue;
        };

        genes_buf.insert((ensembl_id.to_owned(), gene_name.to_owned()));
    }

    Ok(())
}

fn make_enums(
    enum_prefix: &str,
    genes: &HashSet<(String, String)>,
    unavailable_gene_ids: &HashSet<String>,
) -> String {
    let mut ensembl_id_enum_variants = Vec::with_capacity(N_GENES);
    let mut ensembl_id_to_gene_name_match_arms = Vec::with_capacity(N_GENES);

    let ensembl_id_enum_name = Ident::new(&format!("{enum_prefix}EnsemblId"), Span::call_site());

    for (ensembl_id, gene_name) in genes {
        if unavailable_gene_ids.contains(ensembl_id) {
            continue;
        }

        let ensembl_id_variant = Ident::new(&ensembl_id.to_pascal_case(), Span::call_site());

        let gene_name = LitStr::new(gene_name, Span::call_site());

        ensembl_id_to_gene_name_match_arms.push(quote! { Self::#ensembl_id_variant => #gene_name });
        ensembl_id_enum_variants.push(ensembl_id_variant);
    }

    quote! {
        #[derive(Debug, Clone, Copy, ::serde::Deserialize, ::serde::Serialize, PartialEq)]
        #[serde(rename_all = "lowercase")]
        pub enum #ensembl_id_enum_name {
            #(#ensembl_id_enum_variants),*
        }

        impl super::MapToGeneName for #ensembl_id_enum_name {
            #[allow(clippy::match_same_arms, clippy::too_many_lines)]
            fn gene_name(self) -> &'static str {
                match self {
                    #(#ensembl_id_to_gene_name_match_arms),*
                }
            }
        }

    }
    .to_string()
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

fn parse_ensembl_id_and_name_from_gtf_record(record: &RecordBuf) -> Option<(&str, &str)> {
    if record.ty() != "gene" {
        return None;
    }

    let attributes = record.attributes();

    if attributes.get(b"gene_type") != Some(&Value::String(BString::from("protein_coding"))) {
        return None;
    }

    let (Some(Value::String(gene_id)), Some(Value::String(gene_name))) =
        (attributes.get(b"gene_id"), attributes.get(b"gene_name"))
    else {
        unreachable!("'gene_id' and 'gene_name' should be strings");
    };

    Some((
        gene_id
            .to_str()
            .expect("Ensembl ID should be UTF8")
            .split('.')
            .next()
            .unwrap(),
        gene_name.to_str().expect("Gene name should be UTF8"),
    ))
}
