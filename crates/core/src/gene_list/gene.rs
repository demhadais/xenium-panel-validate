mod xenium_prime_human;
mod xenium_prime_mouse;
mod xenium_v1_human;
mod xenium_v1_mouse;

pub trait MapToGeneName {
    fn gene_name(self) -> &'static str;
}
