//! Event types fired by the NJ algorithm and forwarded to wrapper callbacks.

use serde::Serialize;

use crate::alphabet::Alphabet;

/// Log severity level for [`NJEvent::Log`].
#[derive(Clone, Debug, PartialEq, Eq, Serialize, ts_rs::TS)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
pub enum LogLevel {
    Info,
    Warning,
}

/// Events fired by the NJ algorithm and passed to the `on_event` callback.
///
/// Use `#[serde(tag = "type")]` so each variant serialises as a tagged object,
/// e.g. `{ "type": "BootstrapProgress", "completed": 5, "total": 100 }`.
/// This makes it easy to dispatch on `event["type"]` in Python and JavaScript.
#[derive(Clone, Debug, Serialize, ts_rs::TS)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
#[serde(tag = "type")]
pub enum NJEvent {
    /// The MSA passed validation. Fired after input sequences are checked.
    MsaValidated { n_sequences: usize, n_sites: usize },
    /// The alphabet was detected (or confirmed from config).
    AlphabetDetected { alphabet: Alphabet },
    /// Distance matrix computation is about to start.
    ComputingDistances,
    /// Neighbor-Joining algorithm is about to run on the main alignment.
    RunningNJ,
    /// Bootstrap phase is starting. Fired once before the first replicate.
    BootstrapStarted { total: usize },
    /// One bootstrap replicate completed.
    BootstrapProgress { completed: usize, total: usize },
    /// Bootstrap support values are being annotated onto the main tree.
    AnnotatingBootstrap,
    /// A generic log message (warnings, diagnostics, etc.).
    Log { level: LogLevel, message: String },
}
