mod human;
mod mouse;

pub trait MapToGeneName<GN> {
    fn gene_name(self) -> GN;
}

pub trait MapToEnsemblId<EI> {
    fn ensembl_id(self) -> EI;
}

pub trait IsUppercase {
    fn is_uppercase() -> bool;
}

pub trait IsTitlecase {
    fn is_titlecase() -> bool;
}
