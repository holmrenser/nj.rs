//! WebAssembly bindings for the `nj` Neighbor-Joining library.
//!
//! Compiled to `wasm32-unknown-unknown` and bundled via `wasm-pack`. Exposes a
//! single `nj(config)` function to JavaScript. The `config` argument must be a
//! JS object deserialised by `serde-wasm-bindgen` into an [`nj::NJConfig`].
//! Expected shape:
//!
//! ```js
//! {
//!   msa: [{ identifier: "SeqA", sequence: "ACGT" }, ...],
//!   n_bootstrap_samples: 100,
//!   substitution_model: "PDiff",  // "JukesCantor", "Kimura2P", or "Poisson"
//! }
//! ```
//!
//! Returns a Newick string on success, or a JS `Error` string on failure.

use nj::{NJConfig, nj as lib_nj};
use serde_wasm_bindgen::from_value;
use wasm_bindgen::prelude::*;

/// Run Neighbor-Joining and return a Newick string.
///
/// Deserialises `config_json` from a JS object into [`NJConfig`], runs the
/// NJ algorithm, and returns the Newick string. Returns a `JsValue` error
/// string if deserialisation fails or if the NJ algorithm fails (e.g.
/// incompatible model for the detected alphabet, or empty MSA).
#[wasm_bindgen]
pub fn nj(config_json: JsValue) -> Result<String, JsValue> {
    let config: NJConfig = from_value(config_json)
        .map_err(|e| JsValue::from_str(&format!("Invalid NJConfig JSON: {}", e)))?;
    lib_nj(config).map_err(|e| JsValue::from_str(&e.to_string()))
}
