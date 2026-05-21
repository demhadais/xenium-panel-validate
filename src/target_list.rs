use std::path::{Path, PathBuf};

use anyhow::Context;
use serde::Deserialize;

use crate::gene_annotation::GeneAnnotations;

#[derive(Clone, Debug, Deserialize)]
pub struct Target {
    #[serde(alias = "target name", alias = "Target Name")]
    target_name: String,
    #[serde(alias = "ensembl id", alias = "Ensembl ID")]
    ensembl_id: String,
    #[serde(alias = "group", alias = "Group")]
    group: String,
    #[serde(alias = "is backup", alias = "Is Backup")]
    is_backup: bool,
}

pub fn read_target_list_from_path(path: &Path) -> anyhow::Result<Vec<Target>> {
    let context = || format!("failed to read {}", path.to_str().unwrap());

    let csv_reader = csv::Reader::from_path(path).with_context(context)?;

    Ok(csv_reader
        .into_deserialize()
        .collect::<Result<Vec<_>, _>>()
        .with_context(context)?)
}

#[derive(Clone, Debug, thiserror::Error)]
pub enum TargetValidationError {
    #[error("ID {id} not found in genome")]
    TargetIdNotInGenome { id: String },
    #[error(
        "target ID {id} has name {name_in_target_list} in target list, but is called \
         {name_in_genome} in genome"
    )]
    TargetIdNameMismatch {
        id: String,
        name_in_target_list: String,
        name_in_genome: String,
    },
}

pub fn validate_targets_are_in_genome(
    target_list: &[Target],
    gene_annotations: &GeneAnnotations,
) -> Vec<TargetValidationError> {
    let mut errors = Vec::with_capacity(target_list.len());

    for target in target_list {
        if let Err(err) = validate_target_is_in_genome(target, gene_annotations) {
            errors.push(err);
        }
    }

    errors
}

fn validate_target_is_in_genome(
    Target {
        target_name,
        ensembl_id,
        ..
    }: &Target,
    gene_annotations: &GeneAnnotations,
) -> Result<(), TargetValidationError> {
    let Some(target_name_in_genome) = gene_annotations.get(ensembl_id.as_bytes()) else {
        return Err(TargetValidationError::TargetIdNotInGenome {
            id: ensembl_id.to_owned(),
        });
    };

    if target_name.as_bytes() != target_name_in_genome {
        return Err(TargetValidationError::TargetIdNameMismatch {
            id: ensembl_id.to_owned(),
            name_in_target_list: target_name.to_owned(),
            name_in_genome: String::from_utf8(target_name_in_genome.to_owned()).unwrap(),
        });
    }

    Ok(())
}
