mod obs;
mod umi_counts;
mod var;

// pub fn validate_reference_dataset(
//     file: hdf5_metno::File,
//     ReferenceDatasetSpec {
//         cell_annotations_obs_column,
//         ensembl_id_var_column,
//         gene_name_var_column,
//     }: &ReferenceDatasetSpec,
// ) -> Result<ValidatedReferenceDataset, Error> {
//     let mut errors = Vec::with_capacity(16);
//     let mut warnings = Vec::with_capacity(16);

//     let adata = read_reference_dataset(file)?;

//     match check_raw_integer_counts(&adata) {
//         Ok(()) => {}
//         Err(e) => errors.push(e),
//     };

//     match validate_var(&adata, ensembl_id_var_column, gene_name_var_column) {
//         Ok(()) => {}
//         Err(mut e) => errors.append(&mut e),
//     };

//     Ok(ValidatedReferenceDataset {
//         adata,
//         errors,
//         warnings,
//     })
// }

// #[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
// pub struct ReferenceDatasetSpec {
//     cell_annotations_obs_column: String,
//     ensembl_id_var_column: String,
//     gene_name_var_column: String,
// }

// fn validate_counts() {}

// fn check_raw_integer_counts(adata: &AnnData<H5>) -> Result<(), Error> {
//     let x = adata.get_x();

//     let is_integer = x
//         .dtype()
//         .map(|d| d.scalar_type())
//         .flatten()
//         .ok_or(Error::MatrixHasNoDataType)
//         .map(|d| d.is_integer())
//         .unwrap();

//     if !is_integer {
//         return Err(Error::TransformedCounts);
//     }

//     Ok(())
// }

// fn validate_var(
//     adata: &AnnData<H5>,
//     ensembl_id_column: &str,
//     gene_name_column: &str,
// ) -> Result<(), Vec<Error>> {
//     let mut errors = Vec::with_capacity(16);

//     let Ok(var) = adata.get_var().inner().data().map(DataFrame::to_owned) else {
//         errors.push(Error::VarNotFound);
//         return Err(errors);
//     };

//     match check_all_features_have_ensembl_ids(var, ensembl_id_column) {
//         Ok(()) => {}
//         Err(e) => errors.push(e),
//     };

//     if errors.is_empty() {
//         return Ok(());
//     }

//     Err(errors)
// }

// // fn check_all_features_have_ensembl_ids(
// //     var: DataFrame,
// //     ensembl_id_column: &str,
// // ) -> Result<(), Error> {
// //     let ensembl_id_col = col(ensembl_id_column);
// //     let ensembl_id_is_empty_mask = ensembl_id_col
// //         .clone()
// //         .is_null()
// //         .logical_or(ensembl_id_col.str().len_bytes().eq(0));

// //     let ensembl_ids_are_empty = var
// //         .lazy()
// //         .select([ensembl_id_is_empty_mask.alias("ensembl_id_is_empty")])
// //         .collect()
// //         .map_err(|e| Error::EnsemblIdColumnNotFound {
// //             detail: e.to_string(),
// //         })?;

// //     let any_ensembl_id_is_empty = ensembl_ids_are_empty
// //         .column("ensembl_id_is_empty")
// //         .unwrap()
// //         .bool()
// //         .unwrap()
// //         .any();

// //     if any_ensembl_id_is_empty {
// //         return Err(Error::GenesWithoutEnsemblId);
// //     }

// //     Ok(())
// // }

// pub struct ValidatedReferenceDataset {
//     pub adata: AnnData<H5>,
//     pub errors: Vec<Error>,
//     pub warnings: Vec<Warning>,
// }

// #[derive(Clone, Debug, Serialize, PartialEq, Eq)]
// #[serde(tag = "type", rename_all = "snake_case")]
// pub enum Error {
//     InvalidAnnData { reason: String },
//     MatrixHasNoDataType,
//     TransformedCounts,
//     VarNotFound,
//     EnsemblIdColumnNotFound { detail: String },
//     NonstringEnsemblIdColumn { detail: String },
//     GenesWithoutEnsemblId,
// }

// #[derive(Clone, Debug, Serialize, PartialEq, Eq)]
// pub enum Warning {}
