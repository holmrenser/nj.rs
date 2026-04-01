pub use crate::models::SubstitutionModel;
pub use crate::msa::MSA;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, ts_rs::TS, Serialize, Deserialize)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
pub struct SequenceObject {
    pub identifier: String,
    pub sequence: String,
}

impl SequenceObject {
    /// Returns the length of the sequence.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
}

#[derive(Serialize, Deserialize, ts_rs::TS, Clone, Debug)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
pub struct NJConfig {
    pub msa: Vec<SequenceObject>,
    pub n_bootstrap_samples: usize,
    pub substitution_model: SubstitutionModel,
}
