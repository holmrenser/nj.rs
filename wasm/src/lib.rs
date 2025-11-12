use nj::{NJConfig, nj as lib_nj};
use serde_wasm_bindgen::from_value;
use wasm_bindgen::prelude::*;

/// Exposed to JavaScript via wasm-bindgen, this function takes a JSON configuration as a `JsValue`,
/// deserializes it into an `NJConfig`, runs the `nj` algorithm, and returns the result as a `String`.
/// Returns a `JsValue` error if deserialization or computation fails.
#[wasm_bindgen]
pub fn nj(config_json: JsValue) -> Result<String, JsValue> {
    let config: NJConfig = from_value(config_json)
        .map_err(|e| JsValue::from_str(&format!("Invalid NJConfig JSON: {}", e)))?;
    lib_nj(config).map_err(|e| JsValue::from_str(&e.to_string()))
}
