//! Python bindings for the `nj` Neighbor-Joining library.
//!
//! Exposes a single function `nj(config, on_progress=None)` to Python. The
//! `config` argument must be a Python dict (or any dict-like object) that
//! deserialises into an [`nj::NJConfig`] via `serde-pyobject`. Expected shape:
//!
//! ```python
//! {
//!     "msa": [{"identifier": "SeqA", "sequence": "ACGT"}, ...],
//!     "n_bootstrap_samples": 100,
//!     "substitution_model": "PDiff",   # or "JukesCantor", "Kimura2P", "Poisson"
//! }
//! ```
//!
//! `on_progress` is an optional Python callable invoked as
//! `on_progress(completed: int, total: int)` after each bootstrap replicate.
//! It is never called when `n_bootstrap_samples` is 0.
//!
//! Example with tqdm:
//!
//! ```python
//! from tqdm import tqdm
//! from nj_py import nj
//!
//! bar = tqdm(total=100)
//! newick = nj(config, on_progress=lambda current, total: bar.update(1))
//! bar.close()
//! ```
//!
//! Returns a Newick string on success, or raises `ValueError` on invalid
//! config or NJ failure.

#[pyo3::pymodule]
mod _nj_py {
    use ::nj::nj as lib_nj;
    use pyo3::exceptions::PyValueError;
    use pyo3::prelude::*;
    use serde_pyobject::from_pyobject;

    /// Run Neighbor-Joining and return a Newick string.
    ///
    /// `config` must be a dict matching the shape described in the module
    /// documentation. `on_progress`, if provided, is called as
    /// `on_progress(completed, total)` after each bootstrap replicate.
    /// Raises `ValueError` if the config is malformed or if the NJ algorithm
    /// fails (e.g. incompatible model for the detected alphabet).
    #[pyfunction]
    #[pyo3(signature = (py_config, on_progress=None))]
    fn nj(py_config: Bound<PyAny>, on_progress: Option<Py<PyAny>>) -> PyResult<String> {
        let config =
            from_pyobject(py_config).map_err(|e| PyValueError::new_err(e.to_string()))?;

        let callback: Option<Box<dyn Fn(usize, usize)>> = on_progress.map(|cb| {
            Box::new(move |current: usize, total: usize| {
                Python::attach(|py| {
                    cb.call1(py, (current, total)).ok();
                });
            }) as Box<dyn Fn(usize, usize)>
        });

        lib_nj(config, callback).map_err(|e| PyValueError::new_err(e))
    }
}
