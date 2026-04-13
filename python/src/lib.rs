//! Python bindings for the `nj` Neighbor-Joining library.
//!
//! Exposes a single function `nj(config, on_event=None)` to Python. The
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
//! `on_event` is an optional Python callable invoked with a dict describing
//! each algorithm event. The dict always has a `"type"` key identifying the
//! event variant. Bootstrap progress events look like:
//! `{"type": "BootstrapProgress", "completed": 5, "total": 100}`.
//! Stage events include `"MsaValidated"`, `"AlphabetDetected"`,
//! `"ComputingDistances"`, `"RunningNJ"`, `"BootstrapStarted"`,
//! `"AnnotatingBootstrap"`, and `"Log"`.
//!
//! Example with tqdm:
//!
//! ```python
//! from tqdm import tqdm
//! from nj_py import nj
//!
//! bar = None
//! def on_event(event):
//!     global bar
//!     if event["type"] == "BootstrapStarted":
//!         bar = tqdm(total=event["total"])
//!     elif event["type"] == "BootstrapProgress":
//!         bar.update(1)
//!
//! newick = nj(config, on_event=on_event)
//! if bar:
//!     bar.close()
//! ```
//!
//! Returns a Newick string on success, or raises `ValueError` on invalid
//! config or NJ failure.

#[pyo3::pymodule]
mod _nj_py {
    use ::nj::{DistConfig, NJEvent, average_distance as lib_average_distance,
               distance_matrix as lib_distance_matrix, nj as lib_nj};
    use pyo3::exceptions::PyValueError;
    use pyo3::prelude::*;
    use serde_pyobject::{from_pyobject, to_pyobject};

    /// Run Neighbor-Joining and return a Newick string.
    ///
    /// `config` must be a dict matching the shape described in the module
    /// documentation. `on_event`, if provided, is called with a dict for each
    /// algorithm event. See the module documentation for the dict shapes.
    /// Raises `ValueError` if the config is malformed or if the NJ algorithm
    /// fails (e.g. incompatible model for the detected alphabet).
    #[pyfunction]
    #[pyo3(signature = (py_config, on_event=None))]
    fn nj(py_config: Bound<PyAny>, on_event: Option<Py<PyAny>>) -> PyResult<String> {
        let config =
            from_pyobject(py_config).map_err(|e| PyValueError::new_err(e.to_string()))?;

        let callback: Option<Box<dyn Fn(NJEvent)>> = on_event.map(|cb| {
            Box::new(move |event: NJEvent| {
                Python::attach(|py| {
                    if let Ok(obj) = to_pyobject(py, &event) {
                        cb.call1(py, (obj,)).ok();
                    }
                });
            }) as Box<dyn Fn(NJEvent)>
        });

        lib_nj(config, callback).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Compute pairwise distances and return a dict `{"names": [...], "matrix": [[...]]}`.
    ///
    /// `config` must be a dict with shape
    /// `{"msa": [{"identifier": ..., "sequence": ...}, ...], "substitution_model": "PDiff"}`.
    /// Raises `ValueError` if the config is malformed or the model is incompatible.
    #[pyfunction]
    #[pyo3(signature = (py_config))]
    fn distance_matrix(py: Python<'_>, py_config: Bound<PyAny>) -> PyResult<Py<PyAny>> {
        let config: DistConfig =
            from_pyobject(py_config).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let result = lib_distance_matrix(config).map_err(|e| PyValueError::new_err(e.to_string()))?;
        to_pyobject(py, &result)
            .map(|b| b.unbind())
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Compute the mean pairwise distance and return it as a float.
    ///
    /// `config` must be a dict with shape
    /// `{"msa": [{"identifier": ..., "sequence": ...}, ...], "substitution_model": "PDiff"}`.
    /// Raises `ValueError` if the config is malformed or the model is incompatible.
    #[pyfunction]
    #[pyo3(signature = (py_config))]
    fn average_distance(py_config: Bound<PyAny>) -> PyResult<f64> {
        let config: DistConfig =
            from_pyobject(py_config).map_err(|e| PyValueError::new_err(e.to_string()))?;
        lib_average_distance(config).map_err(|e| PyValueError::new_err(e.to_string()))
    }
}
