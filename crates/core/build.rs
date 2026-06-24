use std::{
    collections::HashSet,
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
};

use anyhow::Context;
use bstr::ByteSlice;
use heck::ToPascalCase;
use noodles::gff::feature::{RecordBuf, record_buf::attributes::field::Value};
use proc_macro2::Span;
use quote::quote;
use serde::Deserialize;
use syn::Ident;
use url::Url;

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    println!("cargo::rerun-if-changed=genes.toml");
    let config = include_bytes!("genes.toml");

    let Config { human, mouse } =
        toml::from_slice(config).context("failed to parse config from genes.toml")?;

    // For now, hard-coding these 4 things is no big deal, but in the future, maybe
    // we should refactor this into something reusable
    let human_annotations = || read_gene_annotations(&human.gene_annotations_path);
    let mouse_annotations = || read_gene_annotations(&mouse.gene_annotations_path);

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

    let human_v1_enums = make_enums("XeniumV1Human", human_annotations()?, &human_v1_unavailable)?;
    std::fs::write("src/gene_list/gene/xenium_v1_human.rs", human_v1_enums)?;

    let human_prime_enums = make_enums(
        "XeniumPrimeHuman",
        human_annotations()?,
        &human_prime_unavailable,
    )?;
    std::fs::write(
        "src/gene_list/gene/xenium_prime_human.rs",
        human_prime_enums,
    )?;

    let mouse_v1_enums = make_enums("XeniumV1Mouse", human_annotations()?, &mouse_v1_unavailable)?;
    std::fs::write("src/gene_list/gene/xenium_v1_mouse.rs", mouse_v1_enums)?;

    let mouse_prime_enums = make_enums(
        "XeniumPrimeMouse",
        mouse_annotations()?,
        &mouse_prime_unavailable,
    )?;
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

fn read_gene_annotations(path: &Path) -> anyhow::Result<GeneAnnotationFileReader> {
    let file =
        File::open(path).with_context(|| format!("failed to read {}", path.to_str().unwrap()))?;

    Ok(noodles::gtf::io::Reader::new(std::io::BufReader::new(file)))
}

type GeneAnnotationFileReader = noodles::gtf::io::Reader<BufReader<File>>;

fn make_enums(
    enum_prefix: &str,
    mut reader: GeneAnnotationFileReader,
    unavailable_gene_ids: &HashSet<String>,
) -> anyhow::Result<String> {
    const N_GENES: usize = 30_000;

    let mut ensembl_id_enum_variants = Vec::with_capacity(N_GENES);
    let mut gene_name_enum_variants = Vec::with_capacity(N_GENES);
    let mut ensembl_id_to_gene_name_match_arms = Vec::with_capacity(N_GENES);
    let mut gene_name_to_ensembl_id_match_arms = Vec::with_capacity(N_GENES);

    let ensembl_id_enum_name = Ident::new(&format!("{enum_prefix}EnsemblId"), Span::call_site());
    let gene_name_enum_name = Ident::new(&format!("{enum_prefix}GeneName"), Span::call_site());

    for record in reader.record_bufs() {
        let (ensembl_id, gene_name) = parse_ensembl_id_and_name_from_gtf_record(&record?)?;
        if unavailable_gene_ids.contains(&ensembl_id) {
            continue;
        }

        let ensembl_id_variant = Ident::new(&ensembl_id, Span::call_site());
        let gene_name_variant = Ident::new(&gene_name, Span::call_site());

        ensembl_id_to_gene_name_match_arms
            .push(quote! { Self::#ensembl_id_variant => #gene_name_enum_name::#gene_name_variant });

        gene_name_to_ensembl_id_match_arms.push(
            quote! { Self::#gene_name_variant => #ensembl_id_enum_name::#ensembl_id_variant },
        );

        ensembl_id_enum_variants.push(ensembl_id_variant);
        gene_name_enum_variants.push(gene_name_variant);
    }

    Ok(quote! {
        #[derive(Debug, Clone, Copy, ::serde::Deserialize, ::serde::Serialize, PartialEq)]
        #[serde(rename_all = "lowercase")]
        pub enum #ensembl_id_enum_name {
            #(#ensembl_id_enum_variants),*
        }

        #[derive(Debug, Clone, Copy, ::serde::Deserialize, ::serde::Serialize, PartialEq)]
        #[serde(rename_all = "lowercase")]
        pub enum #gene_name_enum_name {
            #(#gene_name_enum_variants),*
        }

        impl super::MapToGeneName<#gene_name_enum_name> for #ensembl_id_enum_name {
            fn gene_name(self) -> #gene_name_enum_name {
                match self {
                    #(#ensembl_id_to_gene_name_match_arms),*
                }
            }
        }

        impl super::MapToEnsemblId<#ensembl_id_enum_name> for #gene_name_enum_name {
            fn ensembl_id(self) -> #ensembl_id_enum_name {
                match self {
                    #(#gene_name_to_ensembl_id_match_arms),*
                }
            }
        }

    }
    .to_string())
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

fn parse_ensembl_id_and_name_from_gtf_record<'a>(
    record: &'a RecordBuf,
) -> anyhow::Result<(String, String)> {
    let attributes = record.attributes();

    let (Some(Value::String(gene_id)), Some(Value::String(gene_name))) =
        (attributes.get(b"gene_id"), attributes.get(b"gene_name"))
    else {
        unreachable!("'gene_id' and 'gene_name' should be strings");
    };

    Ok((
        gene_id
            .to_str()
            .expect("Ensembl ID should be UTF8")
            .to_pascal_case(),
        gene_name
            .to_str()
            .expect("Ensembl ID should be UTF8")
            .to_pascal_case(),
    ))
}
