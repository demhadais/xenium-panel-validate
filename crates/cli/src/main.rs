use std::fs;

use anyhow::{Context, anyhow};
use camino::{Utf8Path, Utf8PathBuf};
use clap::Parser;
use xenium_panel_validate_cli::targets::{Chemistry, Species, parse_target_list_from_file};

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli {
        Cli::Targets {
            target_path,
            species,
            chemistry,
            field_alias_path,
            field_aliases,
            common:
                CommonOptions {
                    output_format,
                    output,
                },
        } => {
            let parsed_targets = parse_target_list_from_file(
                &target_path,
                species,
                chemistry,
                field_alias_path.as_deref(),
                &field_aliases,
            )?;

            let parsed_targets = match output_format {
                Format::Json => serde_json::to_string(&parsed_targets)?,
            };

            write_report(&parsed_targets, output.as_deref())?;
        }
    }

    Ok(())
}

fn write_report(data: &str, output_path: Option<&Utf8Path>) -> anyhow::Result<()> {
    if let Some(path) = output_path {
        fs::write(path, data).context(format!("failed to write report to {path}"))?;
    } else {
        println!("{data}");
    }

    Ok(())
}

#[derive(clap::Parser)]
enum Cli {
    Targets {
        target_path: Utf8PathBuf,
        #[clap(long, short)]
        species: Species,
        #[clap(long, short)]
        chemistry: Chemistry,
        #[clap(long, short = 'p')]
        field_alias_path: Option<Utf8PathBuf>,
        #[clap(long, short = 'a', value_parser = parse_field_aliases)]
        field_aliases: Vec<(String, String)>,
        #[clap(flatten)]
        common: CommonOptions,
    },
}

fn parse_field_aliases(s: &str) -> anyhow::Result<(String, String)> {
    s.split_once('=')
        .map(|(alias, field)| (alias.to_owned(), field.to_owned()))
        .ok_or_else(|| anyhow!("field aliases must be specified as '<ALIAS>=<FIELD>'"))
}

#[derive(Debug, Clone, PartialEq, Eq, clap::Args)]
struct CommonOptions {
    #[clap(long, short = 'f', default_value_t = Format::Json)]
    output_format: Format,
    #[clap(long, short)]
    output: Option<camino::Utf8PathBuf>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
enum Format {
    Json,
}

impl std::fmt::Display for Format {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Json => "json".fmt(f),
        }
    }
}
