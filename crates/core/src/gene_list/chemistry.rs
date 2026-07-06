#![allow(clippy::unreadable_literal)]
use serde::{Deserialize, Serialize};

mod xenium_prime_human;
mod xenium_prime_mouse;
mod xenium_v1_human;
mod xenium_v1_mouse;

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EnsemblId(&'static str);

impl EnsemblId {
    #[must_use]
    pub fn as_str(&self) -> &'static str {
        self.0
    }
}

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GeneName(&'static str);

impl GeneName {
    #[must_use]
    pub fn as_str(&self) -> &'static str {
        self.0
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct UnvalidatedEnsemblId(String);

impl UnvalidatedEnsemblId {
    #[must_use]
    pub fn new(ensembl_id: String) -> Self {
        Self(ensembl_id)
    }

    #[must_use]
    pub fn to_versionless_uppercase(&self) -> Self {
        Self(self.0.split('.').next().unwrap_or("").to_uppercase())
    }

    #[must_use]
    pub fn is_versionless_and_uppercase(&self) -> bool {
        !self.0.contains('.')
            && self
                .0
                .chars()
                .filter(|c| c.is_alphabetic())
                .all(char::is_uppercase)
    }

    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct UnvalidatedGeneName(String);

impl UnvalidatedGeneName {
    #[must_use]
    pub fn new(gene_name: String) -> Self {
        Self(gene_name)
    }

    #[must_use]
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
pub mod tests {
    use crate::gene_list::chemistry::{
        GeneName, UnvalidatedEnsemblId, xenium_prime_human::XENIUM_PRIME_HUMAN_ENSEMBL_IDS,
        xenium_prime_mouse::XENIUM_PRIME_MOUSE_ENSEMBL_IDS,
        xenium_v1_human::XENIUM_V1_HUMAN_ENSEMBL_IDS, xenium_v1_human_ensembl_id_to_gene_name,
        xenium_v1_mouse::XENIUM_V1_MOUSE_ENSEMBL_IDS,
    };

    pub fn tp53_ensembl_id() -> UnvalidatedEnsemblId {
        UnvalidatedEnsemblId("ENSG00000141510".to_owned())
    }

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
    fn correct_detection_of_versionless_uppercase_ensembl_id() {
        let ensembl_id = tp53_ensembl_id();
        assert!(ensembl_id.is_versionless_and_uppercase());

        let ensembl_id = UnvalidatedEnsemblId(format!("{}.1", ensembl_id.0.to_lowercase()));
        assert!(!ensembl_id.is_versionless_and_uppercase());
    }

    #[test]
    fn canonicalized_ensembl_id_gets_correct_gene_name() {
        let ensembl_id = tp53_ensembl_id().to_versionless_uppercase();

        std::assert_matches!(
            xenium_v1_human_ensembl_id_to_gene_name(&ensembl_id).unwrap(),
            (_, GeneName("TP53"))
        );
    }
}
