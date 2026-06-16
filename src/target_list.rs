use std::path::Path;

use anyhow::Context;
use csv::StringRecord;
use serde::{Deserialize, Serialize};

use crate::gene_annotation::GeneAnnotations;

pub fn read_target_list_from_csv(path: &Path) -> anyhow::Result<Vec<Target>> {
    fn sanitize_row(row: &StringRecord) -> StringRecord {
        StringRecord::from_iter(row.iter().map(|s| s.trim().to_lowercase()))
    }

    let mut input_csv = csv::Reader::from_path(path).context("failed to read target-list")?;
    let headers = input_csv
        .headers()
        .context("failed to read CSV header")
        .map(sanitize_row)?;

    let mut targets = Vec::with_capacity(500);

    for row in input_csv.records() {
        let row = row.map(|row| sanitize_row(&row))?;

        targets.push(
            row.deserialize(Some(&headers))
                .context("failed to deserialize row as target")?,
        );
    }

    Ok(targets)
}

#[derive(Clone, Debug, Deserialize, Serialize, thiserror::Error)]
#[serde(rename_all = "snake_case")]
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

#[derive(Clone, Debug, Deserialize)]
pub struct Target {
    #[serde(alias = "target name")]
    target_name: String,
    #[serde(alias = "ensembl id")]
    ensembl_id: String,
    group: String,
    #[serde(alias = "is backup")]
    is_backup: bool,
}
