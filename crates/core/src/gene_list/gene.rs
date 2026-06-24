mod xenium_prime_human;
mod xenium_prime_mouse;
mod xenium_v1_human;
mod xenium_v1_mouse;

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
