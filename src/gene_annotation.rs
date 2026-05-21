use std::{collections::HashMap, fs::File, io::BufReader, path::Path};

use anyhow::{Context, bail};
use noodles::gff::feature::record_buf::attributes::field::Value;

pub type GeneAnnotationFileReader = noodles::gtf::io::Reader<BufReader<File>>;

pub fn read_gene_annotations(path: &Path) -> anyhow::Result<GeneAnnotationFileReader> {
    let file =
        File::open(path).with_context(|| format!("failed to read {}", path.to_str().unwrap()))?;

    Ok(noodles::gtf::io::Reader::new(std::io::BufReader::new(file)))
}

pub struct GeneAnnotations(HashMap<Vec<u8>, Vec<u8>>);

impl GeneAnnotations {
    pub fn get(&self, gene_id: &[u8]) -> Option<&[u8]> {
        self.0.get(gene_id).map(Vec::as_slice)
    }

    pub fn from_reader(mut reader: GeneAnnotationFileReader) -> anyhow::Result<Self> {
        // This is massive, but realistically shouldn't be that large in terms of memory
        // usage
        let mut gene_annotations = HashMap::with_capacity(2_000_000);

        for (i, record) in reader.record_bufs().enumerate() {
            let record = record.context(format!(
                "failed to read record from line {i} of genes-annotation file"
            ))?;

            let attributes = record.attributes();

            let (Some(Value::String(gene_id)), Some(Value::String(gene_name))) =
                (attributes.get(b"gene_id"), attributes.get(b"gene_name"))
            else {
                bail!(
                    "gene-annotation file does not contain both attributes 'gene_id' and \
                     'gene_name' in line {i} as strings"
                )
            };

            // I love how this library has no way to just own the data, so you have to copy
            // it :)
            gene_annotations.insert(
                gene_id.iter().copied().collect(),
                gene_name.iter().copied().collect(),
            );
        }

        Ok(Self(gene_annotations))
    }
}
