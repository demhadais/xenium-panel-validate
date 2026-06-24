use serde::{Deserialize, Serialize};

use crate::gene_list::gene::{MapToEnsemblId, MapToGeneName};

#[derive(Debug, Clone, Copy, Deserialize, Serialize, PartialEq, Eq, Hash)]
pub enum Grch38EnsemblId {}

impl MapToGeneName<Grch38GeneName> for Grch38EnsemblId {
    fn gene_name(self) -> Grch38GeneName {
        todo!()
    }
}

#[derive(Debug, Clone, Copy, Deserialize, Serialize, PartialEq, Eq, Hash)]
pub enum Grch38GeneName {}

impl MapToEnsemblId<Grch38EnsemblId> for Grch38GeneName {
    fn ensembl_id(self) -> Grch38EnsemblId {
        todo!()
    }
}
