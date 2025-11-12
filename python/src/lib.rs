#[pyo3::pymodule]
mod nj_py {
    use ::nj::nj as lib_nj;
    use pyo3::exceptions::PyValueError;
    use pyo3::prelude::*;
    use serde_pyobject::from_pyobject;

    #[pyfunction]
    fn nj(py_config: Bound<PyAny>) -> PyResult<String> {
        let config = from_pyobject(py_config).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let tree = lib_nj(config).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(tree) // Assuming tree is a String (Newick) or implement ToString
    }
}
