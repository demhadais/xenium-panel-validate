use std::collections::HashMap;

use chemistry::{EnsemblId, GeneName, UnvalidatedEnsemblId, UnvalidatedGeneName};
use csv::StringRecord;
use serde::{Deserialize, Serialize};

pub mod chemistry;

#[allow(clippy::missing_errors_doc)]
pub fn parse_target_list(
    target_list: &str,
    field_aliases: &HashMap<String, String>,
    allowed_genes: &phf::Map<EnsemblId, GeneName>,
) -> csv::Result<ParsedTargetList> {
    let target_list = target_list.trim();
    let mut reader = csv::Reader::from_reader(target_list.as_bytes());

    // We initialize the list of errors from the field-renaming, but it doesn't
    // prevent us from continuing the parsing
    let (fieldnames, error) = rename_fields(reader.headers()?.to_owned(), field_aliases);
    let fieldnames = Some(&fieldnames);

    let mut valid_targets = Vec::with_capacity(500);
    let mut errors = vec![error];

    for (line_number, record) in reader
        .into_records()
        .enumerate()
        .map(|(ln, r)| (Some(ln), r))
    {
        match parse_target_from_record(record?, fieldnames, allowed_genes) {
            Ok(valid_target) => valid_targets.push(valid_target),
            Err(this_records_errors) => errors.push(Error {
                line_number,
                errors: this_records_errors,
            }),
        }
    }

    Ok(ParsedTargetList {
        valid_targets,
        errors,
    })
}

fn rename_fields(
    mut original_fieldnames: StringRecord,
    field_aliases: &HashMap<String, String>,
) -> (StringRecord, Error) {
    original_fieldnames.trim();
    let mut renamed_fields = StringRecord::new();
    let mut errors = Vec::new();

    for original in &original_fieldnames {
        let renamed = field_aliases
            .get(&original.to_lowercase())
            .map_or(original, String::as_str);

        renamed_fields.push_field(renamed);

        if renamed != original {
            errors.push(ErrorInner::RenamedField {
                original_fieldname: original.to_owned(),
                correct_fieldname: renamed.to_owned(),
            });
        }
    }

    (
        renamed_fields,
        Error {
            line_number: None,
            errors,
        },
    )
}

fn parse_target_from_record(
    mut record: StringRecord,
    fieldnames: Option<&StringRecord>,
    allowed_genes: &phf::Map<EnsemblId, GeneName>,
) -> Result<ValidTarget, Vec<ErrorInner>> {
    // Trim the individual fields of the record
    record.trim();

    let Ok(unvalidated_target) = record.deserialize(fieldnames) else {
        // Since every field is an optional string besides the Ensembl ID and gene name,
        // we know that the above deserialization could only fail due to an invalid
        // Ensembl ID or gene name. As such, we report to the user what the invalid
        // values were by deserializing as plain strings
        return Err(vec![ErrorInner::InvalidGene(
            // Unwrapping is fine because extra fields won't cause a failure, nor will missing
            // fields
            record.deserialize(fieldnames).unwrap(),
        )]);
    };

    validate_target(unvalidated_target, allowed_genes)
}

fn validate_target(
    UnvalidatedTarget {
        ensembl_id,
        gene_name,
        group,
        is_backup,
        must_have,
    }: UnvalidatedTarget,
    allowed_genes: &phf::Map<EnsemblId, GeneName>,
) -> Result<ValidTarget, Vec<ErrorInner>> {
    // The number of possible errors in a row is 6 (the same as the number of
    // variants of ErrorInner)
    let mut errors = Vec::with_capacity(8);

    let mut valid_ensembl_id = None;
    let mut valid_gene_name = None;
    let mut valid_is_backup = None;
    let mut valid_must_have = None;

    match validate_ensembl_id_gene_name_pair(ensembl_id.as_ref(), gene_name.as_ref(), allowed_genes)
    {
        Ok((ensembl_id, gene_name)) => {
            valid_ensembl_id = Some(ensembl_id);
            valid_gene_name = Some(gene_name);
        }
        Err(err) => errors.push(err),
    }

    if group.is_none() {
        errors.push(ErrorInner::MissingField("group"));
    }

    match parse_bool_from_str(is_backup.as_deref(), "is_backup") {
        Ok(is_backup) => valid_is_backup = Some(is_backup),
        Err(err) => errors.push(err),
    }

    match parse_bool_from_str(must_have.as_deref(), "must_have") {
        Ok(must_have) => valid_must_have = Some(must_have),
        Err(err) => errors.push(err),
    }

    // Technically, we don't have any compile-time safety that these unwraps are
    // safe. It would be nice to implement that, but I'm not sure how to do that
    // simply at the moment
    if errors.is_empty() {
        Ok(ValidTarget {
            ensembl_id: valid_ensembl_id.unwrap(),
            gene_name: valid_gene_name.unwrap(),
            group: group.unwrap(),
            is_backup: valid_is_backup.unwrap(),
            must_have: valid_must_have.unwrap(),
        })
    } else {
        Err(errors)
    }
}

fn validate_ensembl_id_gene_name_pair(
    ensembl_id: Option<&UnvalidatedEnsemblId>,
    gene_name: Option<&UnvalidatedGeneName>,
    allowed_genes: &phf::Map<EnsemblId, GeneName>,
) -> Result<(EnsemblId, GeneName), ErrorInner> {
    let ensembl_id = ensembl_id.map(UnvalidatedEnsemblId::to_uppercase);

    match (ensembl_id, gene_name) {
        (Some(ensembl_id), maybe_submitted_gene_name) => {
            let (ensembl_id, correct_gene_name) = allowed_genes
                .get_entry(&ensembl_id)
                .map(|(eid, gn)| (*eid, *gn))
                .ok_or_else(|| InvalidGeneError {
                    ensembl_id: Some(ensembl_id),
                    gene_name: maybe_submitted_gene_name.cloned(),
                })?;

            let submitted_gene_name = maybe_submitted_gene_name.ok_or(ErrorInner::NoGeneName {
                ensembl_id,
                probable_gene_name: correct_gene_name,
            })?;

            if *submitted_gene_name != correct_gene_name {
                return Err(ErrorInner::EnsemblIdGeneNameMismatch {
                    ensembl_id,
                    submitted_gene_name: submitted_gene_name.to_owned(),
                    correct_gene_name,
                });
            }

            Ok((ensembl_id, correct_gene_name))
        }

        (None, Some(gene_name)) => Err(ErrorInner::NoEnsemblId {
            gene_name: gene_name.to_owned(),
        }),
        (None, None) => Err(ErrorInner::MissingGene),
    }
}

fn parse_bool_from_str(s: Option<&str>, fieldname: &'static str) -> Result<bool, ErrorInner> {
    let Some(s) = s else {
        return Err(ErrorInner::MissingField(fieldname));
    };

    s.parse().map_err(|_| ErrorInner::ParseBool {
        value: s.to_owned(),
    })
}

#[derive(Clone, Debug, Serialize)]
pub struct ParsedTargetList {
    valid_targets: Vec<ValidTarget>,
    errors: Vec<Error>,
}

#[derive(Clone, Debug, Serialize)]
struct ValidTarget {
    ensembl_id: EnsemblId,
    gene_name: GeneName,
    group: String,
    is_backup: bool,
    must_have: bool,
}

#[derive(Clone, Debug, Deserialize)]
struct UnvalidatedTarget {
    ensembl_id: Option<UnvalidatedEnsemblId>,
    gene_name: Option<UnvalidatedGeneName>,
    group: Option<String>,
    is_backup: Option<String>,
    must_have: Option<String>,
}

#[derive(Clone, Debug, Serialize, PartialEq)]
struct Error {
    line_number: Option<usize>,
    errors: Vec<ErrorInner>,
}

#[derive(Clone, Debug, Serialize, PartialEq)]
enum ErrorInner {
    MissingGene,
    MissingField(&'static str),
    ParseBool {
        value: String,
    },
    NoEnsemblId {
        gene_name: UnvalidatedGeneName,
    },
    NoGeneName {
        ensembl_id: EnsemblId,
        probable_gene_name: GeneName,
    },
    RenamedField {
        original_fieldname: String,
        correct_fieldname: String,
    },
    EnsemblIdGeneNameMismatch {
        ensembl_id: EnsemblId,
        submitted_gene_name: UnvalidatedGeneName,
        correct_gene_name: GeneName,
    },
    // TODO: Unknown or Unrecognized?
    InvalidGene(InvalidGeneError),
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
struct InvalidGeneError {
    ensembl_id: Option<UnvalidatedEnsemblId>,
    gene_name: Option<UnvalidatedGeneName>,
}

impl From<InvalidGeneError> for ErrorInner {
    fn from(err: InvalidGeneError) -> Self {
        Self::InvalidGene(err)
    }
}
