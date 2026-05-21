use clap::Parser;
use xenium_prepaneldesign_validator::{TargetListValidationSettings, validate_target_list};

#[derive(Parser)]
enum Cli {
    ValidateTargetList(TargetListValidationSettings),
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    let errors = match cli {
        Cli::ValidateTargetList(settings) => validate_target_list(&settings)?,
    };

    if errors.len() != 0 {
        dbg!(&errors[0]);
    }

    Ok(())
}
