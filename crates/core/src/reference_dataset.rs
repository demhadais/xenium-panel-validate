// Commment for commit
use std::path::{Path, PathBuf};

use anndata::{
    AnnData, AnnDataOp, ArrayElemOp, Backend,
    backend::DataType,
    container::{Inner, InnerDataFrameElem},
};
use anndata_hdf5::{H5, H5File};
use polars::prelude::*;
use serde::{Deserialize, Serialize};

mod nonempty_anndata;

pub fn validate_reference_dataset(
    file: H5File,
    ReferenceDatasetSpec {
        cell_annotations_obs_column,
        ensembl_id_var_column,
        gene_name_var_column,
    }: &ReferenceDatasetSpec,
) -> Result<ValidatedReferenceDataset, Error> {
    let mut errors = Vec::with_capacity(16);
    let mut warnings = Vec::with_capacity(16);

    let adata = read_reference_dataset(file)?;

    match check_raw_integer_counts(&adata) {
        Ok(()) => {}
        Err(e) => errors.push(e),
    };

    match validate_var(&adata, ensembl_id_var_column, gene_name_var_column) {
        Ok(()) => {}
        Err(mut e) => errors.append(&mut e),
    };

    Ok(ValidatedReferenceDataset {
        adata,
        errors,
        warnings,
    })
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub struct ReferenceDatasetSpec {
    cell_annotations_obs_column: String,
    ensembl_id_var_column: String,
    gene_name_var_column: String,
}

fn read_reference_dataset(file: H5File) -> Result<AnnData<H5>, Error> {
    match AnnData::open(file) {
        Ok(a) => Ok(a),
        Err(e) => Err(Error::InvalidAnnData {
            reason: e.to_string(),
        }),
    }
}

fn check_raw_integer_counts(adata: &AnnData<H5>) -> Result<(), Error> {
    let x = adata.get_x();

    let is_integer = x
        .dtype()
        .map(|d| d.scalar_type())
        .flatten()
        .ok_or(Error::MatrixHasNoDataType)
        .map(|d| d.is_integer())
        .unwrap();

    if !is_integer {
        return Err(Error::TransformedCounts);
    }

    Ok(())
}

fn validate_var(
    adata: &AnnData<H5>,
    ensembl_id_column: &str,
    gene_name_column: &str,
) -> Result<(), Vec<Error>> {
    let mut errors = Vec::with_capacity(16);

    let Ok(var) = adata.get_var().inner().data().map(DataFrame::to_owned) else {
        errors.push(Error::VarNotFound);
        return Err(errors);
    };

    match check_all_features_have_ensembl_ids(var, ensembl_id_column) {
        Ok(()) => {}
        Err(e) => errors.push(e),
    };

    if errors.is_empty() {
        return Ok(());
    }

    Err(errors)
}

fn check_all_features_have_ensembl_ids(
    var: DataFrame,
    ensembl_id_column: &str,
) -> Result<(), Error> {
    let ensembl_id_col = col(ensembl_id_column);
    let ensembl_id_is_empty_mask = ensembl_id_col
        .clone()
        .is_null()
        .logical_or(ensembl_id_col.str().len_bytes().eq(0));

    let ensembl_ids_are_empty = var
        .lazy()
        .select([ensembl_id_is_empty_mask.alias("ensembl_id_is_empty")])
        .collect()
        .map_err(|e| Error::EnsemblIdColumnNotFound {
            detail: e.to_string(),
        })?;

    let any_ensembl_id_is_empty = ensembl_ids_are_empty
        .column("ensembl_id_is_empty")
        .unwrap()
        .bool()
        .unwrap()
        .any();

    if any_ensembl_id_is_empty {
        return Err(Error::GenesWithoutEnsemblId);
    }

    Ok(())
}

pub struct ValidatedReferenceDataset {
    pub adata: AnnData<H5>,
    pub errors: Vec<Error>,
    pub warnings: Vec<Warning>,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Error {
    InvalidAnnData { reason: String },
    MatrixHasNoDataType,
    TransformedCounts,
    VarNotFound,
    EnsemblIdColumnNotFound { detail: String },
    NonstringEnsemblIdColumn { detail: String },
    GenesWithoutEnsemblId,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub enum Warning {}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use anndata::Backend;
    use anndata_hdf5::{H5, H5File};

    use crate::reference_dataset::{
        ReferenceDatasetSpec, read_reference_dataset, validate_reference_dataset,
    };

    #[test]
    fn valid_anndata() {
        let path = Path::new(
            "test-data/SOD1_G93A_mouse_spinal_cord_P112_specimen_1_SOD1_G93A_mouse_spinal_cord_P112_specimen_1_sample_filtered_feature_bc_matrix_small.h5ad",
        );

        let file = H5::open(path).unwrap();

        let report = validate_reference_dataset(
            file,
            &ReferenceDatasetSpec {
                cell_annotations_obs_column: "".to_owned(),
                ensembl_id_var_column: "gene_id".to_owned(),
                gene_name_var_column: "gene_name".to_owned(),
            },
        )
        .unwrap();

        assert_eq!(report.errors, []);
    }
}
