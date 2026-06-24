pub mod xenium_prime_human;
pub mod xenium_prime_mouse;
pub mod xenium_v1_human;
pub mod xenium_v1_mouse;

pub trait MapToGeneName {
    fn gene_name(self) -> &'static str;
}
