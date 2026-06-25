use std::collections::HashMap;

use chemistry::{EnsemblId, GeneName, UnvalidatedEnsemblId, UnvalidatedGeneName};
use csv::StringRecord;
use serde::{Deserialize, Serialize};

pub mod chemistry;

#[allow(clippy::missing_errors_doc)]
pub fn parse_target_list(
    target_list: &str,
    field_aliases: &HashMap<String, String>,
    allowed_genes: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)> + Copy,
) -> csv::Result<ParsedTargetList> {
    let target_list = target_list.trim();
    let mut reader = csv::Reader::from_reader(target_list.as_bytes());

    // We initialize the list of errors from the field-renaming, but it doesn't
    // prevent us from continuing the parsing
    let (fieldnames, error) = rename_fields(reader.headers()?.to_owned(), field_aliases);
    let fieldnames = Some(&fieldnames);
    let mut errors = vec![error];

    let mut valid_targets = Vec::with_capacity(500);

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
    allowed_genes: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
) -> Result<ValidTarget, Vec<ErrorInner>> {
    // Trim the individual fields of the record
    record.trim();

    // Unwrapping is fine because extra fields won't cause a failure, nor will missing fields
    let unvalidated_target = record.deserialize(fieldnames).unwrap();

    validate_target(unvalidated_target, allowed_genes)
}

fn validate_target(
    UnvalidatedTarget {
        gene,
        group,
        is_backup,
        must_have,
    }: UnvalidatedTarget,
    allowed_genes: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
) -> Result<ValidTarget, Vec<ErrorInner>> {
    // The number of possible errors in a row is 6 (the same as the number of
    // variants of ErrorInner)
    let mut errors = Vec::with_capacity(8);

    let mut valid_ensembl_id = None;
    let mut valid_gene_name = None;
    let mut valid_is_backup = None;
    let mut valid_must_have = None;

    match validate_ensembl_id_gene_name_pair(&gene, allowed_genes) {
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
    unvalidated_gene: &UnvalidatedGene,
    allowed_genes: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
) -> Result<(EnsemblId, GeneName), ErrorInner> {
    let UnvalidatedGene {
        ensembl_id,
        gene_name,
    } = unvalidated_gene;

    let ensembl_id = ensembl_id.as_ref().map(UnvalidatedEnsemblId::to_uppercase);

    match (ensembl_id, gene_name.as_ref()) {
        (Some(ensembl_id), maybe_submitted_gene_name) => {
            let (ensembl_id, correct_gene_name) =
                allowed_genes(&ensembl_id).ok_or_else(|| unvalidated_gene.to_owned())?;

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

#[derive(Clone, Debug, Deserialize, PartialEq)]
struct UnvalidatedTarget {
    #[serde(flatten)]
    gene: UnvalidatedGene,
    group: Option<String>,
    is_backup: Option<String>,
    must_have: Option<String>,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
struct UnvalidatedGene {
    ensembl_id: Option<UnvalidatedEnsemblId>,
    gene_name: Option<UnvalidatedGeneName>,
}

#[derive(Clone, Debug, Serialize, PartialEq)]
struct Error {
    line_number: Option<usize>,
    errors: Vec<ErrorInner>,
}

#[derive(Clone, Debug, Serialize, PartialEq)]
#[serde(rename_all = "snake_case", tag = "type")]
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
    InvalidGene(UnvalidatedGene),
}

impl From<UnvalidatedGene> for ErrorInner {
    fn from(err: UnvalidatedGene) -> Self {
        Self::InvalidGene(err)
    }
}

#[cfg(test)]
mod tests {
    use crate::gene_list::{
        Error, ErrorInner, UnvalidatedGene, UnvalidatedTarget,
        chemistry::{
            UnvalidatedEnsemblId, UnvalidatedGeneName, xenium_v1_human_ensembl_id_to_gene_name,
        },
        rename_fields, validate_ensembl_id_gene_name_pair,
    };

    #[test]
    fn renaming_fields() {
        let original_fieldnames = ["field1", "field2"].iter().collect();
        let field_aliases = [("field1", "field_1")]
            .map(|(k, v)| (k.to_owned(), v.to_owned()))
            .into_iter()
            .collect();

        let (renamed_fields, error) = rename_fields(original_fieldnames, &field_aliases);

        assert_eq!(
            renamed_fields,
            ["field_1", "field2"][..],
            "failed to rename fields"
        );

        assert_eq!(
            error,
            Error {
                line_number: None,
                errors: vec![ErrorInner::RenamedField {
                    original_fieldname: "field1".to_owned(),
                    correct_fieldname: "field_1".to_owned()
                }]
            },
            "failed to construct field-renaming error"
        );
    }

    #[test]
    fn unvalidated_target_deserializes_with_invalid_fields() {
        let data = b"field1,field2,ensembl_id\nvalue1,value2,id";
        let mut reader = csv::Reader::from_reader(&data[..]);

        let deserialized: Vec<UnvalidatedTarget> =
            reader.deserialize().collect::<Result<_, _>>().unwrap();

        assert_eq!(
            deserialized,
            vec![UnvalidatedTarget {
                gene: UnvalidatedGene {
                    ensembl_id: Some(UnvalidatedEnsemblId("id".to_owned())),
                    gene_name: None
                },
                group: None,
                is_backup: None,
                must_have: None
            }]
        )
    }

    fn tp53_ensembl_id() -> UnvalidatedEnsemblId {
        UnvalidatedEnsemblId("ENSG00000141510".to_owned())
    }

    #[test]
    fn valid_gene() {
        validate_ensembl_id_gene_name_pair(
            &UnvalidatedGene {
                ensembl_id: Some(tp53_ensembl_id()),
                gene_name: Some(UnvalidatedGeneName("TP53".to_owned())),
            },
            xenium_v1_human_ensembl_id_to_gene_name,
        )
        .unwrap();
    }

    #[test]
    fn ensembl_id_gene_name_mismatch() {
        let ensembl_id = tp53_ensembl_id();
        let gene_name = UnvalidatedGeneName(String::new());

        let err = validate_ensembl_id_gene_name_pair(
            &UnvalidatedGene {
                ensembl_id: Some(ensembl_id.clone()),
                gene_name: Some(gene_name.clone()),
            },
            xenium_v1_human_ensembl_id_to_gene_name,
        )
        .unwrap_err();

        let (correct_ensembl_id, correct_gene_name) =
            xenium_v1_human_ensembl_id_to_gene_name(&ensembl_id).unwrap();

        assert_eq!(
            err,
            ErrorInner::EnsemblIdGeneNameMismatch {
                ensembl_id: correct_ensembl_id,
                submitted_gene_name: gene_name,
                correct_gene_name,
            },
            "failed to create Ensembl ID-gene name mismatch error"
        );
    }
}
