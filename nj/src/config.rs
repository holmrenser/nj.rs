pub use crate::msa::MSA;
use serde::{Deserialize, Serialize};
use ts_rs::TS;

#[derive(Clone, Debug, TS, Serialize, Deserialize)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
pub struct FastaSequence {
    pub identifier: String,
    pub sequence: String,
}

impl FastaSequence {
    /// Returns the length of the sequence.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
}

#[derive(Serialize, Deserialize, TS, Clone, Debug)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
pub struct NJConfig {
    pub msa: Vec<FastaSequence>,
    pub show_internal: bool,
    pub n_bootstrap_samples: usize,
}
