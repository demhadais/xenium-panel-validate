use std::collections::{HashMap, HashSet};

use chemistry::{EnsemblId, GeneName, UnvalidatedEnsemblId, UnvalidatedGeneName};
use csv::StringRecord;
use serde::{Deserialize, Serialize};

pub mod chemistry;

#[allow(clippy::missing_errors_doc)]
pub fn parse_target_list(
    target_list: &str,
    field_aliases: &HashMap<&str, &str>,
    ensembl_id_to_gene: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)> + Copy,
) -> csv::Result<ParsedTargetList> {
    const N_GENES: usize = 500;
    let target_list = target_list.trim();
    let mut reader = csv::Reader::from_reader(target_list.as_bytes());

    // We initialize the list of errors from the field-renaming, but it doesn't
    // prevent us from continuing the parsing
    let (fieldnames, error) = rename_fields(reader.headers()?.to_owned(), field_aliases);
    let fieldnames = Some(&fieldnames);
    let mut errors = error.map(|e| vec![e]).unwrap_or_default();

    let mut valid_targets = Vec::with_capacity(N_GENES);
    let mut genes = HashSet::with_capacity(N_GENES);

    for (line_number, record) in reader
        .into_records()
        .enumerate()
        .map(|(ln, r)| (Some(ln), r))
    {
        let (submitted_target, result) =
            parse_target_from_record(record?, fieldnames, ensembl_id_to_gene);
        let submitted_target = Some(submitted_target);

        match result {
            Ok(vt) => {
                let is_new = genes.insert(vt.gene);

                if is_new {
                    valid_targets.push(vt);
                } else {
                    errors.push(Error {
                        line_number,
                        submitted_target,
                        errors: vec![ErrorInner::DuplicateGene],
                    });
                }
            }
            Err(errs) => errors.push(Error {
                line_number,
                submitted_target,
                errors: errs,
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
    field_aliases: &HashMap<&str, &str>,
) -> (StringRecord, Option<Error>) {
    original_fieldnames.trim();
    let mut renamed_fields = StringRecord::new();
    let mut errors = Vec::new();

    for original in &original_fieldnames {
        let renamed = field_aliases
            .get(original.to_lowercase().as_str())
            .unwrap_or(&original);

        renamed_fields.push_field(renamed);

        if *renamed != original {
            errors.push(ErrorInner::RenamedField {
                original_fieldname: original.to_owned(),
                correct_fieldname: (*renamed).to_owned(),
            });
        }
    }

    (
        renamed_fields,
        (!errors.is_empty()).then_some(Error {
            line_number: None,
            submitted_target: None,
            errors,
        }),
    )
}

#[allow(clippy::result_large_err)]
fn parse_target_from_record(
    mut record: StringRecord,
    fieldnames: Option<&StringRecord>,
    ensembl_id_to_gene_name: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
) -> (UnvalidatedTarget, Result<ValidTarget, Vec<ErrorInner>>) {
    // Trim the individual fields of the record
    record.trim();

    // Unwrapping is fine because extra fields won't cause a failure, nor will
    // missing fields
    let unvalidated_target = record.deserialize(fieldnames).unwrap();
    let validation_result = validate_target(&unvalidated_target, ensembl_id_to_gene_name);

    (unvalidated_target, validation_result)
}

fn validate_target(
    UnvalidatedTarget {
        gene,
        group,
        is_backup,
        must_have,
    }: &UnvalidatedTarget,
    ensembl_id_to_gene_name: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
) -> Result<ValidTarget, Vec<ErrorInner>> {
    let mut errors = Vec::with_capacity(4);

    if group.is_none() {
        errors.push(ErrorInner::MissingField("group"));
    }

    let is_backup = match parse_bool_from_str(is_backup.as_deref(), "is_backup") {
        Ok(is_backup) => Some(is_backup),
        Err(err) => {
            errors.push(err);
            None
        }
    };

    let must_have = match parse_bool_from_str(must_have.as_deref(), "must_have") {
        Ok(must_have) => Some(must_have),
        Err(err) => {
            errors.push(err);
            None
        }
    };

    if let (Some(is_backup), Some(must_have)) = (is_backup, must_have)
        && is_backup
        && must_have
    {
        errors.push(ErrorInner::BackupAndMustHave);
    }

    let valid_gene = match validate_ensembl_id_gene_name_pair(gene, ensembl_id_to_gene_name) {
        Ok(vg) => Some(vg),
        Err(err) => {
            errors.push(err);
            None
        }
    };

    match (valid_gene, group, is_backup, must_have) {
        (Some(valid_gene), Some(group), Some(is_backup), Some(must_have)) => Ok(ValidTarget {
            gene: valid_gene,
            group: group.to_lowercase(),
            is_backup,
            must_have,
        }),
        _ => Err(errors),
    }
}

fn validate_ensembl_id_gene_name_pair(
    unvalidated_gene: &UnvalidatedGene,
    ensembl_id_to_gene_name: impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)>,
) -> Result<ValidGene, ErrorInner> {
    let UnvalidatedGene {
        ensembl_id,
        gene_name: submitted_gene_name,
    } = unvalidated_gene;

    let Some(ensembl_id) = ensembl_id else {
        return Err(ErrorInner::NoEnsemblId);
    };

    let map_valid_gene = |(ensembl_id, gene_name)| ValidGene {
        ensembl_id,
        gene_name,
    };

    let valid_gene = if ensembl_id.is_versionless_and_uppercase() {
        ensembl_id_to_gene_name(ensembl_id)
            .map(map_valid_gene)
            .ok_or(ErrorInner::GeneNotFound)?
    } else {
        let maybe_valid_gene =
            ensembl_id_to_gene_name(&ensembl_id.to_versionless_uppercase()).map(map_valid_gene);

        return Err(ErrorInner::VersionedOrLowercaseEnsemblId {
            correct_gene: maybe_valid_gene,
        });
    };

    let Some(submitted_gene_name) = submitted_gene_name else {
        return Err(ErrorInner::NoGeneName {
            probable_gene_name: valid_gene.gene_name,
        });
    };

    if *submitted_gene_name == valid_gene.gene_name {
        Ok(valid_gene)
    } else {
        Err(ErrorInner::EnsemblIdGeneNameMismatch {
            correct_gene_name: valid_gene.gene_name,
        })
    }
}

fn parse_bool_from_str(s: Option<&str>, fieldname: &'static str) -> Result<bool, ErrorInner> {
    let Some(s) = s else {
        return Err(ErrorInner::MissingField(fieldname));
    };

    let s = s.to_lowercase();

    s.parse().map_err(|_| ErrorInner::ParseBool { value: s })
}

#[derive(Clone, Debug, Serialize)]
pub struct ParsedTargetList {
    pub valid_targets: Vec<ValidTarget>,
    pub errors: Vec<Error>,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct ValidTarget {
    #[serde(flatten)]
    gene: ValidGene,
    group: String,
    is_backup: bool,
    must_have: bool,
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, Hash)]
pub struct ValidGene {
    ensembl_id: EnsemblId,
    gene_name: GeneName,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
struct UnvalidatedTarget {
    #[serde(flatten)]
    gene: UnvalidatedGene,
    group: Option<String>,
    is_backup: Option<String>,
    must_have: Option<String>,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
struct UnvalidatedGene {
    ensembl_id: Option<UnvalidatedEnsemblId>,
    gene_name: Option<UnvalidatedGeneName>,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct Error {
    line_number: Option<usize>,
    submitted_target: Option<UnvalidatedTarget>,
    errors: Vec<ErrorInner>,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum ErrorInner {
    MissingField(&'static str),
    ParseBool {
        value: String,
    },
    VersionedOrLowercaseEnsemblId {
        correct_gene: Option<ValidGene>,
    },
    NoEnsemblId,
    NoGeneName {
        probable_gene_name: GeneName,
    },
    RenamedField {
        original_fieldname: String,
        correct_fieldname: String,
    },
    EnsemblIdGeneNameMismatch {
        correct_gene_name: GeneName,
    },
    BackupAndMustHave,
    GeneNotFound,
    DuplicateGene,
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::gene_list::{
        Error, ErrorInner, ParsedTargetList, UnvalidatedGene, UnvalidatedTarget, ValidGene,
        ValidTarget,
        chemistry::{
            UnvalidatedEnsemblId, UnvalidatedGeneName, tests::tp53_ensembl_id,
            xenium_v1_human_ensembl_id_to_gene_name,
        },
        parse_target_list, rename_fields, validate_ensembl_id_gene_name_pair,
    };

    #[test]
    fn renaming_fields() {
        let original_fieldnames = ["field1", "field2"].iter().collect();
        let field_aliases = [("field1", "field_1")].into_iter().collect();

        let (renamed_fields, error) = rename_fields(original_fieldnames, &field_aliases);

        assert_eq!(
            renamed_fields,
            ["field_1", "field2"][..],
            "failed to rename fields"
        );

        assert_eq!(
            error,
            Some(Error {
                line_number: None,
                submitted_target: None,
                errors: vec![ErrorInner::RenamedField {
                    original_fieldname: "field1".to_owned(),
                    correct_fieldname: "field_1".to_owned()
                }]
            }),
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
            [UnvalidatedTarget {
                gene: UnvalidatedGene {
                    ensembl_id: Some(UnvalidatedEnsemblId::new("id".to_owned())),
                    gene_name: None
                },
                group: None,
                is_backup: None,
                must_have: None
            }]
        )
    }

    #[test]
    fn valid_gene() {
        let _ = validate_ensembl_id_gene_name_pair(
            &UnvalidatedGene {
                ensembl_id: Some(tp53_ensembl_id()),
                gene_name: Some(UnvalidatedGeneName::new("TP53".to_owned())),
            },
            xenium_v1_human_ensembl_id_to_gene_name,
        )
        .unwrap();
    }

    #[test]
    fn ensembl_id_gene_name_mismatch() {
        let ensembl_id = tp53_ensembl_id();
        let gene_name = UnvalidatedGeneName::new(String::new());

        let err = validate_ensembl_id_gene_name_pair(
            &UnvalidatedGene {
                ensembl_id: Some(ensembl_id.clone()),
                gene_name: Some(gene_name.clone()),
            },
            xenium_v1_human_ensembl_id_to_gene_name,
        )
        .unwrap_err();

        let (_correct_ensembl_id, correct_gene_name) =
            xenium_v1_human_ensembl_id_to_gene_name(&ensembl_id).unwrap();

        assert_eq!(
            err,
            ErrorInner::EnsemblIdGeneNameMismatch { correct_gene_name },
            "failed to create Ensembl ID-gene name mismatch error"
        );
    }

    #[test]
    fn duplicate_targets() {
        let ensembl_id = tp53_ensembl_id();
        let ensembl_id_str = ensembl_id.as_str();

        // Two rows with the same Ensembl ID/gene-name pair but differing other fields
        let gene_list = format!(
            "ensembl_id,gene_name,group,is_backup,must_have\n{ensembl_id_str},TP53,group0,false,\
             true\n{ensembl_id_str},TP53,group1,true,false"
        );

        let ParsedTargetList {
            valid_targets,
            errors,
        } = parse_target_list(
            &gene_list,
            &HashMap::new(),
            xenium_v1_human_ensembl_id_to_gene_name,
        )
        .unwrap();

        let (correct_eid, correct_gn) =
            xenium_v1_human_ensembl_id_to_gene_name(&ensembl_id).unwrap();

        assert_eq!(
            valid_targets,
            [ValidTarget {
                gene: ValidGene {
                    ensembl_id: correct_eid,
                    gene_name: correct_gn
                },
                group: "group0".to_owned(),
                is_backup: false,
                must_have: true
            }]
        );

        assert_eq!(
            errors.as_array::<1>().unwrap()[0].errors,
            [ErrorInner::DuplicateGene],
            "did not find exactly 1 error: {:?}",
            errors
        );
    }
}
