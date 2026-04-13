//! WebAssembly bindings for the `nj` Neighbor-Joining library.
//!
//! Compiled to `wasm32-unknown-unknown` and bundled via `wasm-pack`. Exposes a
//! single `nj(config, onEvent?)` function to JavaScript. The `config`
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
//! `onEvent` is an optional JS function called with a tagged event object for
//! each algorithm stage. The object always has a `type` field identifying the
//! event. Bootstrap progress events look like:
//! `{ type: 'BootstrapProgress', completed: 5, total: 100 }`.
//!
//! ```js
//! import { nj } from '@holmrenser/nj';
//!
//! const newick = nj(config, (event) => {
//!   if (event.type === 'BootstrapProgress') {
//!     progressBar.value = event.completed / event.total * 100;
//!   }
//! });
//! ```
//!
//! Returns a Newick string on success, or throws a JS error string on failure.

use js_sys::Function;
use nj::{DistConfig, NJConfig, NJEvent, average_distance as lib_average_distance,
         distance_matrix as lib_distance_matrix, nj as lib_nj};
use serde_wasm_bindgen::{from_value, to_value};
use wasm_bindgen::prelude::*;

/// Run Neighbor-Joining and return a Newick string.
///
/// Deserialises `config_json` from a JS object into [`NJConfig`], runs the
/// NJ algorithm, and returns the Newick string. `on_event`, if provided,
/// is called with a tagged event object at each algorithm stage. See the
/// module documentation for the event object shapes.
/// Returns a `JsValue` error string if deserialisation fails or if the NJ
/// algorithm fails (e.g. incompatible model for the detected alphabet, or
/// empty MSA).
#[wasm_bindgen]
pub fn nj(config_json: JsValue, on_event: Option<Function>) -> Result<String, JsValue> {
    let config: NJConfig = from_value(config_json)
        .map_err(|e| JsValue::from_str(&format!("Invalid NJConfig: {}", e)))?;

    let callback: Option<Box<dyn Fn(NJEvent)>> = on_event.map(|f| {
        Box::new(move |event: NJEvent| {
            if let Ok(js_val) = to_value(&event) {
                let args = js_sys::Array::of1(&js_val);
                f.apply(&JsValue::NULL, &args).ok();
            }
        }) as Box<dyn Fn(NJEvent)>
    });

    lib_nj(config, callback).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Compute pairwise distances and return a `{ names, matrix }` JS object.
///
/// `config_json` must be a JS object with shape:
/// `{ msa: [{ identifier, sequence }, ...], substitution_model: 'PDiff' }`.
/// Throws a JS error string if the config is invalid or the model is incompatible.
#[wasm_bindgen]
pub fn distance_matrix(config_json: JsValue) -> Result<JsValue, JsValue> {
    let config: DistConfig = from_value(config_json)
        .map_err(|e| JsValue::from_str(&format!("Invalid DistConfig: {}", e)))?;
    let result = lib_distance_matrix(config).map_err(|e| JsValue::from_str(&e.to_string()))?;
    to_value(&result).map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e)))
}

/// Compute the mean pairwise distance and return it as a JS number.
///
/// `config_json` must be a JS object with shape:
/// `{ msa: [{ identifier, sequence }, ...], substitution_model: 'PDiff' }`.
/// Throws a JS error string if the config is invalid or the model is incompatible.
#[wasm_bindgen]
pub fn average_distance(config_json: JsValue) -> Result<f64, JsValue> {
    let config: DistConfig = from_value(config_json)
        .map_err(|e| JsValue::from_str(&format!("Invalid DistConfig: {}", e)))?;
    lib_average_distance(config).map_err(|e| JsValue::from_str(&e.to_string()))
}
