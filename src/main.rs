use std::path::{Path, PathBuf};

use anyhow::Context;
use clap::Parser;
use serde::Serialize;
use xenium_prepaneldesign_validator::{TargetListValidationSettings, validate_target_list};

#[derive(Parser)]
enum Cli {
    ValidateTargets {
        #[clap(flatten)]
        settings: TargetListValidationSettings,
        #[clap(short, long)]
        errors_path: PathBuf,
    },
}

fn write_errors(path: &Path, errors: &[impl Serialize]) -> anyhow::Result<()> {
    std::fs::write(path, serde_json::to_string(errors)?).context(format!(
        "failed to write errors to {}",
        path.to_str().unwrap()
    ))?;

    Ok(())
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    let (errors_path, errors) = match cli {
        Cli::ValidateTargets {
            settings,
            errors_path,
        } => (errors_path, validate_target_list(&settings)?),
    };

    write_errors(&errors_path, &errors)?;

    Ok(())
}
