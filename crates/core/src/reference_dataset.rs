// Commment for commit
use std::path::{Path, PathBuf};

use anndata::{AnnData, ArrayElemOp, Backend, backend::DataType};
use anndata_hdf5::{H5, H5File};
use serde::{Deserialize, Serialize};

mod anndata_helpers;

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

    match has_raw_integer_counts(&adata) {
        Ok(true) => {}
        Ok(false) => errors.push(Error::TransformedCounts),
        Err(e) => errors.push(e),
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

fn has_raw_integer_counts(adata: &AnnData<H5>) -> Result<bool, Error> {
    let x = adata.get_x();

    if x.is_none() {
        return Ok(false);
    }

    x.dtype()
        .map(|d| d.scalar_type())
        .flatten()
        .ok_or(Error::MatrixHasNoDataType)
        .map(|d| d.is_integer())
}

fn all_genes_have_ensembl_ids(adata: &AnnData<H5>) -> Result<bool, Error> {
    let var = adata.get_var();
    if var.is_none() {
        return Err(Error::NoVarAnnotations);
    }

    Ok(todo!())
}

pub struct ValidatedReferenceDataset {
    pub adata: AnnData<H5>,
    pub errors: Vec<Error>,
    pub warnings: Vec<Warning>,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub enum Error {
    InvalidAnnData { reason: String },
    MatrixHasNoDataType,
    TransformedCounts,
    NoVarAnnotations,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub enum Warning {}
