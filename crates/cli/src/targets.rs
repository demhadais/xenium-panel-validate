use std::{collections::HashMap, fs};

use anyhow::Context;
use camino::Utf8Path;
use xenium_panel_validate::gene_list::{
    ParsedTargetList,
    chemistry::{
        xenium_prime_human_ensembl_id_to_gene_name, xenium_prime_mouse_ensembl_id_to_gene_name,
        xenium_v1_human_ensembl_id_to_gene_name, xenium_v1_mouse_ensembl_id_to_gene_name,
    },
    parse_target_list,
};

pub fn parse_target_list_from_file(
    target_path: &Utf8Path,
    species: Species,
    chemistry: Chemistry,
    field_alias_path: Option<&Utf8Path>,
    field_aliases: &[(String, String)],
) -> anyhow::Result<ParsedTargetList> {
    let target_list = fs::read_to_string(target_path)
        .context(format!("failed to read target-list from {target_path}"))?;

    let field_alias_file_contents =
        field_alias_path
            .map(fs::read)
            .transpose()
            .with_context(|| {
                if let Some(field_alias_path) = field_alias_path {
                    format!("failed to read field aliases from {field_alias_path}")
                } else {
                    "failed to read field aliases".to_owned()
                }
            })?;

    let field_aliases = read_field_aliases(field_alias_file_contents.as_deref(), &field_aliases)
        .with_context(|| {
            if let Some(path) = field_alias_path {
                format!("failed to read field aliases from {path}",)
            } else {
                format!("failed to construct field aliases")
            }
        })?;

    // This repetition sucks, but the only way to fix it is to make a closure that
    // takes a Box<impl Fn(&UnvalidatedEnsemblId) -> Option<(EnsemblId, GeneName)> +
    // Copy>, which is not worth it for just 4 repetitions
    let result = match (species, chemistry) {
        (Species::HomoSapiens, Chemistry::V1) => parse_target_list(
            &target_list,
            &field_aliases,
            xenium_v1_human_ensembl_id_to_gene_name,
        ),
        (Species::HomoSapiens, Chemistry::Prime) => parse_target_list(
            &target_list,
            &field_aliases,
            xenium_prime_human_ensembl_id_to_gene_name,
        ),
        (Species::MusMusculus, Chemistry::V1) => parse_target_list(
            &target_list,
            &field_aliases,
            xenium_v1_mouse_ensembl_id_to_gene_name,
        ),
        (Species::MusMusculus, Chemistry::Prime) => parse_target_list(
            &target_list,
            &field_aliases,
            xenium_prime_mouse_ensembl_id_to_gene_name,
        ),
    };

    Ok(result?)
}

fn read_field_aliases<'a>(
    field_alias_file_contents: Option<&'a [u8]>,
    field_aliases: &'a [(String, String)],
) -> anyhow::Result<HashMap<&'a str, &'a str>> {
    let mut field_aliases: HashMap<_, _> = field_aliases
        .iter()
        .map(|(s1, s2)| (s1.as_str(), s2.as_str()))
        .collect();

    let Some(aliases_from_file) = field_alias_file_contents else {
        return Ok(field_aliases);
    };

    let aliases_from_file: HashMap<_, _> = toml::from_slice(&aliases_from_file)?;

    for (alias, field) in aliases_from_file {
        // We want field-aliases from the CLI to take precedence
        if !field_aliases.contains_key(alias) {
            field_aliases.insert(alias, field);
        }
    }

    Ok(field_aliases)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum Species {
    HomoSapiens,
    MusMusculus,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum Chemistry {
    V1,
    Prime,
}

#[cfg(test)]
mod tests {
    use crate::targets::read_field_aliases;

    #[test]
    fn field_aliases_are_combined_correctly() {
        let field_aliases = ["alias1", "field1", "alias2", "field2"];

        let field_aliases: Vec<(String, String)> = field_aliases
            .chunks(2)
            .map(|alias_field| (alias_field[0].to_owned(), alias_field[1].to_owned()))
            .collect();

        let field_aliases =
            read_field_aliases(Some(br#"alias1 = "field2""#), &field_aliases).unwrap();

        assert_eq!(
            field_aliases,
            [("alias1", "field1"), ("alias2", "field2")]
                .into_iter()
                .collect()
        );
    }
}
