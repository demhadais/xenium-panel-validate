use phf::{PhfEq, PhfHash};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EnsemblId(&'static str);

impl PhfHash for EnsemblId {
    fn phf_hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.0.phf_hash(state);
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GeneName(&'static str);

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct UnvalidatedEnsemblId(String);

impl UnvalidatedEnsemblId {
    pub fn to_uppercase(&self) -> Self {
        Self(self.0.to_uppercase())
    }

    pub fn get(&self) -> &str {
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
pub struct UnvalidatedGeneName(String);

impl PartialEq<GeneName> for UnvalidatedGeneName {
    fn eq(&self, other: &GeneName) -> bool {
        self.0 == other.0
    }
}

impl PartialEq<UnvalidatedEnsemblId> for GeneName {
    fn eq(&self, other: &UnvalidatedEnsemblId) -> bool {
        self.0 == other.0
    }
}
