//! WebAssembly bindings for the `nj` Neighbor-Joining library.
//!
//! Compiled to `wasm32-unknown-unknown` and bundled via `wasm-pack`. Exposes a
//! single `nj(config, onProgress?)` function to JavaScript. The `config`
//! argument must be a JS object deserialised by `serde-wasm-bindgen` into an
//! [`nj::NJConfig`]. Expected shape:
//!
//! ```js
//! {
//!   msa: [{ identifier: 'SeqA', sequence: 'ACGT' }, ...],
//!   n_bootstrap_samples: 100,
//!   substitution_model: 'PDiff',  // 'JukesCantor', 'Kimura2P', or 'Poisson'
//! }
//! ```
//!
//! `onProgress` is an optional JS function called as `onProgress(completed, total)`
//! after each bootstrap replicate. It is never called when
//! `n_bootstrap_samples` is 0.
//!
//! ```js
//! import { nj } from '@holmrenser/nj';
//!
//! const newick = nj(config, (current, total) => {
//!   progressBar.value = current / total * 100;
//! });
//! ```
//!
//! Returns a Newick string on success, or throws a JS error string on failure.

use js_sys::Function;
use nj::{NJConfig, nj as lib_nj};
use serde_wasm_bindgen::from_value;
use wasm_bindgen::prelude::*;

/// Run Neighbor-Joining and return a Newick string.
///
/// Deserialises `config_json` from a JS object into [`NJConfig`], runs the
/// NJ algorithm, and returns the Newick string. `on_progress`, if provided,
/// is called as `on_progress(completed, total)` after each bootstrap replicate.
/// Returns a `JsValue` error string if deserialisation fails or if the NJ
/// algorithm fails (e.g. incompatible model for the detected alphabet, or
/// empty MSA).
#[wasm_bindgen]
pub fn nj(config_json: JsValue, on_progress: Option<Function>) -> Result<String, JsValue> {
    let config: NJConfig = from_value(config_json)
        .map_err(|e| JsValue::from_str(&format!("Invalid NJConfig: {}", e)))?;

    let callback: Option<Box<dyn Fn(usize, usize)>> = on_progress.map(|f| {
        Box::new(move |current: usize, total: usize| {
            let args = js_sys::Array::of2(
                &JsValue::from(current as u32),
                &JsValue::from(total as u32),
            );
            f.apply(&JsValue::NULL, &args).ok();
        }) as Box<dyn Fn(usize, usize)>
    });

    lib_nj(config, callback).map_err(|e| JsValue::from_str(&e))
}
