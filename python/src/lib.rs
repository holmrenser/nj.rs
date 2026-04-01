//! Python bindings for the `nj` Neighbor-Joining library.
//!
//! Exposes a single function `nj(config)` to Python. The `config` argument
//! must be a Python dict (or any dict-like object) that deserialises into an
//! [`nj::NJConfig`] via `serde-pyobject`. Expected shape:
//!
//! ```python
//! {
//!     "msa": [{"identifier": "SeqA", "sequence": "ACGT"}, ...],
//!     "n_bootstrap_samples": 100,
//!     "substitution_model": "PDiff",   # or "JukesCantor", "Kimura2P", "Poisson"
//! }
//! ```
//!
//! Returns a Newick string on success, or raises `ValueError` on invalid
//! config or NJ failure.

#[pyo3::pymodule]
mod nj_py {
    use ::nj::nj as lib_nj;
    use pyo3::exceptions::PyValueError;
    use pyo3::prelude::*;
    use serde_pyobject::from_pyobject;

    /// Run Neighbor-Joining and return a Newick string.
    ///
    /// `config` must be a dict matching the shape described in the module
    /// documentation. Raises `ValueError` if the config is malformed or if the
    /// NJ algorithm fails (e.g. incompatible model for the detected alphabet).
    #[pyfunction]
    fn nj(py_config: Bound<PyAny>) -> PyResult<String> {
        let config = from_pyobject(py_config).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let tree = lib_nj(config).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(tree)
    }
}
