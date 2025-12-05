use serde::{Deserialize, Serialize};
use ts_rs::TS;

#[derive(Clone, Debug, TS, Serialize, Deserialize)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
pub struct FastaSequence {
    pub identifier: String,
    pub sequence: String,
}

pub type MSA = Vec<FastaSequence>;

#[derive(Serialize, Deserialize, TS, Clone, Debug)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
pub struct NJConfig {
    pub msa: MSA,
    pub hide_internal: bool,
}
