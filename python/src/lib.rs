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
//!     # optional: include extra fields in the result
//!     "include_distance_matrix": True,
//!     "include_average_distance": True,
//! }
//! ```
//!
//! Returns a dict with a `"newick"` key always present. When
//! `"include_distance_matrix"` is `True`, the dict also contains a
//! `"distance_matrix"` key (`{"names": [...], "matrix": [[...]]}`). When
//! `"include_average_distance"` is `True`, the dict also contains an
//! `"average_distance"` key (a float). Both flags may be combined.
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
//! result = nj(config, on_event=on_event)
//! print(result["newick"])
//! if bar:
//!     bar.close()
//! ```
//!
//! Raises `ValueError` on invalid config or NJ failure.

#[pyo3::pymodule]
mod _nj_py {
    use ::nj::{DistConfig, NJError, NJEvent, average_distance as lib_average_distance,
               distance_matrix as lib_distance_matrix, nj as lib_nj};
    use pyo3::exceptions::{PyRuntimeError, PyValueError};
    use pyo3::prelude::*;
    use serde_pyobject::{from_pyobject, to_pyobject};

    /// Maps a library [`NJError`] to an appropriate Python exception.
    ///
    /// Usage/validation errors (bad input, incompatible model) become
    /// `ValueError`; unexpected internal failures become `RuntimeError`. The
    /// `NJError::code()` is prefixed so callers can still branch on the exact
    /// variant if they need finer granularity than the exception class.
    fn map_nj_error(e: NJError) -> PyErr {
        let msg = format!("[{}] {}", e.code(), e);
        match e {
            NJError::EmptyMsa
            | NJError::EmptySequence
            | NJError::SequenceLengthMismatch { .. }
            | NJError::DuplicateIdentifier { .. }
            | NJError::IncompatibleModel { .. }
            | NJError::ParseError(_) => PyValueError::new_err(msg),
            NJError::AlgorithmFailure(_) | NJError::RngError(_) => PyRuntimeError::new_err(msg),
        }
    }

    /// Run Neighbor-Joining and return a result dict.
    ///
    /// `config` must be a dict matching the shape described in the module
    /// documentation. `on_event`, if provided, is called with a dict for each
    /// algorithm event. See the module documentation for the dict shapes.
    /// Raises `ValueError` if the config is malformed or if the NJ algorithm
    /// fails (e.g. incompatible model for the detected alphabet).
    ///
    /// The returned dict always has a `"newick"` key. It also has a
    /// `"distance_matrix"` key when `include_distance_matrix` is `True`, and
    /// an `"average_distance"` key when `include_average_distance` is `True`.
    #[pyfunction]
    #[pyo3(signature = (py_config, on_event=None))]
    fn nj(py: Python<'_>, py_config: Bound<PyAny>, on_event: Option<Py<PyAny>>) -> PyResult<Py<PyAny>> {
        let config =
            from_pyobject(py_config).map_err(|e| PyValueError::new_err(e.to_string()))?;

        let callback: Option<Box<dyn Fn(NJEvent)>> = on_event.map(|cb| {
            Box::new(move |event: NJEvent| {
                Python::attach(|py| match to_pyobject(py, &event) {
                    // Surface a failing user callback's traceback to stderr
                    // instead of silently swallowing it.
                    Ok(obj) => {
                        if let Err(err) = cb.call1(py, (obj,)) {
                            err.print(py);
                        }
                    }
                    Err(e) => {
                        eprintln!("nj_py: failed to serialize event for callback: {e}");
                    }
                });
            }) as Box<dyn Fn(NJEvent)>
        });

        let result = lib_nj(config, callback).map_err(map_nj_error)?;
        to_pyobject(py, &result)
            .map(|b| b.unbind())
            .map_err(|e| PyValueError::new_err(e.to_string()))
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
        let result = lib_distance_matrix(config).map_err(map_nj_error)?;
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
        lib_average_distance(config).map_err(map_nj_error)
    }
}
