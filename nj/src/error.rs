//! Error type for the `nj` library.

use crate::alphabet::Alphabet;
use crate::models::SubstitutionModel;
use std::fmt;

/// Errors that can occur during Neighbor-Joining analysis.
#[derive(Debug)]
pub enum NJError {
    /// The input MSA contains no sequences.
    EmptyMsa,
    /// A sequence in the MSA has zero length.
    EmptySequence,
    /// Sequences in the MSA have different lengths.
    SequenceLengthMismatch {
        expected: usize,
        got: usize,
        identifier: String,
    },
    /// Two or more sequences share the same identifier. Duplicate names would
    /// corrupt bootstrap clade counting (clades are keyed by name → index).
    DuplicateIdentifier { identifier: String },
    /// The substitution model is incompatible with the detected or specified alphabet.
    IncompatibleModel {
        model: SubstitutionModel,
        alphabet: Alphabet,
    },
    /// The gamma rate-heterogeneity shape parameter is not finite and positive.
    InvalidGammaShape { value: f64 },
    /// The proportion of invariant sites is not finite and in `[0, 1)`.
    InvalidPInvar { value: f64 },
    /// An internal NJ algorithm failure (should be unreachable for valid input).
    AlgorithmFailure(String),
    /// Failed to collect entropy from the OS PRNG (bootstrap only).
    RngError(String),
    /// Failed to parse a FASTA-formatted string.
    ParseError(String),
}

impl NJError {
    /// Returns a stable, machine-readable code identifying the error variant.
    ///
    /// Intended for language bindings that want to expose the error kind
    /// programmatically (e.g. as a JavaScript `Error.name`) rather than relying
    /// on the human-readable [`Display`](std::fmt::Display) message.
    pub fn code(&self) -> &'static str {
        match self {
            NJError::EmptyMsa => "EmptyMsa",
            NJError::EmptySequence => "EmptySequence",
            NJError::SequenceLengthMismatch { .. } => "SequenceLengthMismatch",
            NJError::DuplicateIdentifier { .. } => "DuplicateIdentifier",
            NJError::IncompatibleModel { .. } => "IncompatibleModel",
            NJError::InvalidGammaShape { .. } => "InvalidGammaShape",
            NJError::InvalidPInvar { .. } => "InvalidPInvar",
            NJError::AlgorithmFailure(_) => "AlgorithmFailure",
            NJError::RngError(_) => "RngError",
            NJError::ParseError(_) => "ParseError",
        }
    }
}

impl fmt::Display for NJError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NJError::EmptyMsa => write!(f, "Input MSA is empty"),
            NJError::EmptySequence => write!(f, "Sequences must not be empty"),
            NJError::SequenceLengthMismatch { expected, got, identifier } => write!(
                f,
                "All sequences must have the same length. Expected {expected}, got {got} for '{identifier}'"
            ),
            NJError::DuplicateIdentifier { identifier } => {
                write!(f, "Duplicate sequence identifier: '{identifier}'")
            }
            NJError::IncompatibleModel { model, alphabet } => write!(
                f,
                "Substitution model {model:?} is incompatible with {alphabet:?} alphabet"
            ),
            NJError::InvalidGammaShape { value } => write!(
                f,
                "Gamma shape parameter must be a finite positive number, got {value}"
            ),
            NJError::InvalidPInvar { value } => write!(
                f,
                "Proportion of invariant sites must be a finite number in [0, 1), got {value}"
            ),
            NJError::AlgorithmFailure(msg) => write!(f, "NJ algorithm failure: {msg}"),
            NJError::RngError(msg) => write!(f, "RNG error: {msg}"),
            NJError::ParseError(msg) => write!(f, "FASTA parse error: {msg}"),
        }
    }
}

impl std::error::Error for NJError {}
