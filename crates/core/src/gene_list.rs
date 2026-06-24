use std::collections::HashMap;

use csv::StringRecord;
use serde::{Deserialize, Serialize, de::DeserializeOwned};

use crate::gene_list::gene::MapToGeneName;

pub mod gene;

#[allow(clippy::missing_errors_doc)]
pub fn parse_target_list<EI>(
    target_list: &str,
    field_aliases: &HashMap<String, String>,
) -> csv::Result<ParsedTargetList<EI>>
where
    EI: MapToGeneName + DeserializeOwned + Copy,
{
    // Trim and lowercase the whole target-list to avoid whitespace and casing
    // errors
    let target_list = target_list.trim().to_lowercase();
    let mut reader = csv::Reader::from_reader(target_list.as_bytes());

    // We initialize the list of errors from the field-renaming, but it doesn't
    // prevent us from continuing the parsing
    let (fieldnames, error) = rename_fields(reader.headers()?, field_aliases);
    let fieldnames = Some(&fieldnames);

    let mut valid_targets = Vec::with_capacity(500);
    let mut errors = vec![error];

    for (line_number, record) in reader
        .into_records()
        .enumerate()
        .map(|(ln, r)| (Some(ln), r))
    {
        match parse_target_from_record(record?, fieldnames) {
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

fn rename_fields<EI>(
    original_fieldnames: &StringRecord,
    field_aliases: &HashMap<String, String>,
) -> (StringRecord, Error<EI>) {
    let mut renamed_fields = StringRecord::new();
    let mut errors = Vec::new();

    for original in original_fieldnames {
        let renamed = field_aliases.get(original).map_or(original, String::as_str);

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

fn parse_target_from_record<EI>(
    mut record: StringRecord,
    fieldnames: Option<&StringRecord>,
) -> Result<ValidTarget<EI>, Vec<ErrorInner<EI>>>
where
    EI: MapToGeneName + DeserializeOwned + Copy,
{
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

    validate_target(unvalidated_target)
}

fn validate_target<EI>(
    UnvalidatedTarget {
        ensembl_id,
        gene_name,
        group,
        is_backup,
        must_have,
    }: UnvalidatedTarget<EI>,
) -> Result<ValidTarget<EI>, Vec<ErrorInner<EI>>>
where
    EI: MapToGeneName + Copy,
{
    // The number of possible errors in a row is 6 (the same as the number of
    // variants of ErrorInner)
    let mut errors = Vec::with_capacity(8);

    let mut valid_ensembl_id = None;
    let mut valid_gene_name = None;
    let mut valid_is_backup = None;
    let mut valid_must_have = None;

    match validate_ensembl_id_gene_name_pair(ensembl_id, gene_name.as_deref()) {
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

fn validate_ensembl_id_gene_name_pair<EI>(
    ensembl_id: Option<EI>,
    gene_name: Option<&str>,
) -> Result<(EI, &'static str), ErrorInner<EI>>
where
    EI: MapToGeneName + Copy,
{
    match (ensembl_id, gene_name) {
        (Some(ensembl_id), Some(submitted_gene_name)) => {
            let correct_gene_name = ensembl_id.gene_name();

            if correct_gene_name == submitted_gene_name {
                Ok((ensembl_id, correct_gene_name))
            } else {
                Err(ErrorInner::EnsemblIdGeneNameMismatch {
                    ensembl_id,
                    submitted_gene_name: submitted_gene_name.to_owned(),
                    correct_gene_name,
                })
            }
        }
        (Some(ensembl_id), None) => Err(ErrorInner::NoGeneName {
            ensembl_id,
            probable_gene_name: ensembl_id.gene_name(),
        }),
        (None, Some(gene_name)) => Err(ErrorInner::NoEnsemblId {
            gene_name: gene_name.to_owned(),
        }),
        (None, None) => Err(ErrorInner::MissingGene),
    }
}

fn parse_bool_from_str<EI>(
    s: Option<&str>,
    fieldname: &'static str,
) -> Result<bool, ErrorInner<EI>> {
    let Some(s) = s else {
        return Err(ErrorInner::MissingField(fieldname));
    };

    s.parse().map_err(|_| ErrorInner::ParseBool {
        value: s.to_owned(),
    })
}

#[derive(Clone, Debug, Serialize)]
pub struct ParsedTargetList<EI> {
    valid_targets: Vec<ValidTarget<EI>>,
    errors: Vec<Error<EI>>,
}

#[derive(Clone, Debug, Serialize)]
struct ValidTarget<EI> {
    ensembl_id: EI,
    gene_name: &'static str,
    group: String,
    is_backup: bool,
    must_have: bool,
}

#[derive(Clone, Debug, Deserialize)]
struct UnvalidatedTarget<EI> {
    ensembl_id: Option<EI>,
    gene_name: Option<String>,
    group: Option<String>,
    is_backup: Option<String>,
    must_have: Option<String>,
}

#[derive(Clone, Debug, Serialize, PartialEq)]
struct Error<EI> {
    line_number: Option<usize>,
    errors: Vec<ErrorInner<EI>>,
}

#[derive(Clone, Debug, Serialize, PartialEq)]
enum ErrorInner<EI> {
    MissingGene,
    MissingField(&'static str),
    ParseBool {
        value: String,
    },
    NoEnsemblId {
        gene_name: String,
    },
    NoGeneName {
        ensembl_id: EI,
        probable_gene_name: &'static str,
    },
    RenamedField {
        original_fieldname: String,
        correct_fieldname: String,
    },
    EnsemblIdGeneNameMismatch {
        ensembl_id: EI,
        submitted_gene_name: String,
        correct_gene_name: &'static str,
    },
    // TODO: Unknown or Unrecognized?
    InvalidGene(InvalidGeneError),
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
struct InvalidGeneError {
    ensembl_id: Option<String>,
    gene_name: Option<String>,
}
