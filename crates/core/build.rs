use std::{fs::File, io::BufReader, path::Path};

use anyhow::Context;
use noodles::gff::feature::record_buf::attributes::field::Value;
use quote::quote;

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    println!("cargo::rerun-if-changed=genes.toml");
    Ok(())
}

fn read_gene_annotations(path: &Path) -> anyhow::Result<GeneAnnotationFileReader> {
    let file =
        File::open(path).with_context(|| format!("failed to read {}", path.to_str().unwrap()))?;

    Ok(noodles::gtf::io::Reader::new(std::io::BufReader::new(file)))
}

type GeneAnnotationFileReader = noodles::gtf::io::Reader<BufReader<File>>;

fn make_enums(reader: GeneAnnotationFileReader) -> proc_macro2::TokenStream {
    quote! {}
}
