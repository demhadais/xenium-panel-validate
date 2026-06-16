use std::path::PathBuf;

use anyhow::Context;
use serde::{Deserialize, Serialize};

use crate::{
    gene_annotation::{GeneAnnotations, read_gene_annotations},
    target_list::{
        TargetValidationError, read_target_list_from_csv, validate_targets_are_in_genome,
    },
};

mod gene_annotation;
mod target_list;

#[derive(Debug, Serialize, Deserialize, clap::Args)]
pub struct TargetListValidationSettings {
    #[clap(short, long)]
    target_list_path: PathBuf,
    #[clap(short, long)]
    gene_annotations_path: PathBuf,
}

#[allow(clippy::missing_errors_doc)]
pub fn validate_target_list(
    TargetListValidationSettings {
        target_list_path,
        gene_annotations_path,
    }: &TargetListValidationSettings,
) -> anyhow::Result<Vec<TargetValidationError>> {
    let target_list =
        read_target_list_from_csv(target_list_path).context("failed to read target list")?;
    let gene_annotations_reader =
        read_gene_annotations(gene_annotations_path).context("failed to read gene annotations")?;

    let genes = GeneAnnotations::from_reader(gene_annotations_reader)
        .context("failed to parse gene annotations")?;

    let errors = validate_targets_are_in_genome(&target_list, &genes);

    Ok(errors)
}
