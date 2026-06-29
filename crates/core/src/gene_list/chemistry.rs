#![allow(clippy::unreadable_literal)]
use phf::{PhfEq, PhfHash};
use serde::{Deserialize, Serialize};

mod xenium_prime_human;
mod xenium_prime_mouse;
mod xenium_v1_human;
mod xenium_v1_mouse;

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EnsemblId(&'static str);

impl EnsemblId {
    pub fn as_str(&self) -> &'static str {
        self.0
    }
}

impl PhfHash for EnsemblId {
    fn phf_hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.0.phf_hash(state);
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GeneName(&'static str);

impl GeneName {
    pub fn as_str(&self) -> &'static str {
        self.0
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct UnvalidatedEnsemblId(pub String);

impl UnvalidatedEnsemblId {
    #[must_use]
    pub fn to_versionless_uppercased(&self) -> Self {
        Self(self.0.split('.').next().unwrap().to_uppercase())
    }

    #[must_use]
    pub fn is_versionless_and_uppercase(&self) -> bool {
        !self.0.contains('.')
            && self
                .0
                .chars()
                .filter(|c| c.is_alphanumeric())
                .all(char::is_uppercase)
    }

    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl PhfHash for UnvalidatedEnsemblId {
    fn phf_hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.0.phf_hash(state);
    }
}

impl PhfEq<UnvalidatedEnsemblId> for EnsemblId {
    fn phf_eq(&self, other: &UnvalidatedEnsemblId) -> bool {
        self.0.phf_eq(&other.0)
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct UnvalidatedGeneName(pub String);

impl UnvalidatedGeneName {
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl PartialEq<GeneName> for UnvalidatedGeneName {
    fn eq(&self, other: &GeneName) -> bool {
        self.0 == other.0
    }
}

#[must_use]
pub fn xenium_v1_human_ensembl_id_to_gene_name(
    ensembl_id: &UnvalidatedEnsemblId,
) -> Option<(EnsemblId, GeneName)> {
    ensembl_id_to_gene_name(ensembl_id, &xenium_v1_human::XENIUM_V1_HUMAN_ENSEMBL_IDS)
}

#[must_use]
pub fn xenium_prime_human_ensembl_id_to_gene_name(
    ensembl_id: &UnvalidatedEnsemblId,
) -> Option<(EnsemblId, GeneName)> {
    ensembl_id_to_gene_name(
        ensembl_id,
        &xenium_prime_human::XENIUM_PRIME_HUMAN_ENSEMBL_IDS,
    )
}

#[must_use]
pub fn xenium_v1_mouse_ensembl_id_to_gene_name(
    ensembl_id: &UnvalidatedEnsemblId,
) -> Option<(EnsemblId, GeneName)> {
    ensembl_id_to_gene_name(ensembl_id, &xenium_v1_mouse::XENIUM_V1_MOUSE_ENSEMBL_IDS)
}

#[must_use]
pub fn xenium_prime_mouse_ensembl_id_to_gene_name(
    ensembl_id: &UnvalidatedEnsemblId,
) -> Option<(EnsemblId, GeneName)> {
    ensembl_id_to_gene_name(
        ensembl_id,
        &xenium_prime_mouse::XENIUM_PRIME_MOUSE_ENSEMBL_IDS,
    )
}

fn ensembl_id_to_gene_name(
    ensembl_id: &UnvalidatedEnsemblId,
    map: &phf::Map<&'static str, &'static str>,
) -> Option<(EnsemblId, GeneName)> {
    map.get_entry(&ensembl_id.0)
        .map(|(eid, gn)| (EnsemblId(eid), GeneName(gn)))
}

#[cfg(test)]
mod tests {
    use crate::gene_list::chemistry::{
        GeneName, UnvalidatedEnsemblId, xenium_prime_human::XENIUM_PRIME_HUMAN_ENSEMBL_IDS,
        xenium_prime_mouse::XENIUM_PRIME_MOUSE_ENSEMBL_IDS,
        xenium_v1_human::XENIUM_V1_HUMAN_ENSEMBL_IDS, xenium_v1_human_ensembl_id_to_gene_name,
        xenium_v1_mouse::XENIUM_V1_MOUSE_ENSEMBL_IDS,
    };

    #[test]
    fn unavailable_genes_are_not_in_map() {
        let unavailable_genes = [
            ("ENSG00000273816", &XENIUM_V1_HUMAN_ENSEMBL_IDS),
            ("ENSG00000249966", &XENIUM_PRIME_HUMAN_ENSEMBL_IDS),
            ("ENSMUSG00000117061", &XENIUM_V1_MOUSE_ENSEMBL_IDS),
            ("ENSMUSG00000094028", &XENIUM_PRIME_MOUSE_ENSEMBL_IDS),
        ];

        for (ensembl_id, map) in unavailable_genes {
            std::assert_matches!(map.get(ensembl_id), None);
        }
    }

    #[test]
    fn canonicalized_ensembl_id_gets_correct_gene_name() {
        let ensembl_id =
            UnvalidatedEnsemblId("eNSG00000141510.1".to_string()).to_versionless_uppercased();

        std::assert_matches!(
            xenium_v1_human_ensembl_id_to_gene_name(&ensembl_id).unwrap(),
            (_, GeneName("TP53"))
        );
    }
}
