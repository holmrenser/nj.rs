//! Public configuration types shared between the core library, Python bindings,
//! and WASM bindings.
//!
//! Both [`SequenceObject`] and [`NJConfig`] implement `Serialize`/`Deserialize`
//! so they can be round-tripped through JSON (WASM), Python dicts (`serde-pyobject`),
//! and TypeScript via the generated `lib_types.ts` type definitions.

pub use crate::alphabet::Alphabet;
pub use crate::models::SubstitutionModel;
pub use crate::msa::MSA;
use serde::{Deserialize, Serialize};

/// A single sequence in a multiple sequence alignment.
///
/// Both fields are plain `String`s so that the struct can be constructed
/// from any source (FASTA parser, JSON, Python dict) without encoding overhead.
/// The actual byte-level encoding into typed symbols is deferred to
/// [`crate::msa::MSA::push`].
#[derive(Clone, Debug, ts_rs::TS, Serialize, Deserialize)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
pub struct SequenceObject {
    /// Sequence name, e.g. the FASTA header without the leading `>`.
    pub identifier: String,
    /// Raw sequence string, e.g. `"ACGT"` or `"ARND"`. Gap characters (`-`)
    /// are preserved and handled by the substitution models.
    pub sequence: String,
}

impl SequenceObject {
    /// Returns the length of the sequence in bytes (not symbols).
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
}

/// Configuration for distance-only computation (no NJ, no bootstrap).
///
/// Pass a `DistConfig` to [`crate::distance_matrix`] or [`crate::average_distance`].
/// The alphabet (DNA vs. protein) is auto-detected from the sequences unless
/// [`alphabet`](DistConfig::alphabet) is explicitly set;
/// [`substitution_model`](DistConfig::substitution_model) must be compatible with
/// the alphabet or an error is returned.
#[derive(Serialize, Deserialize, ts_rs::TS, Clone, Debug)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
pub struct DistConfig {
    /// The aligned sequences. All sequences must have the same length.
    pub msa: Vec<SequenceObject>,
    /// Substitution model used to compute pairwise distances.
    pub substitution_model: SubstitutionModel,
    /// Override automatic alphabet detection. When `None` (the default), the
    /// alphabet is inferred from the sequences. Set to `Some(Alphabet::DNA)` or
    /// `Some(Alphabet::Protein)` to force a specific alphabet.
    #[serde(default)]
    pub alphabet: Option<Alphabet>,
    /// Maximum number of threads to use for parallel distance computation.
    ///
    /// Only effective when the `parallel` Cargo feature is enabled.
    /// When `None` (the default), Rayon uses all available hardware threads.
    #[serde(default)]
    pub num_threads: Option<usize>,
}

/// Full configuration for a single Neighbor-Joining run.
///
/// Pass an `NJConfig` to [`crate::nj`] to run the algorithm and receive a
/// Newick string. The alphabet (DNA vs. protein) is auto-detected from the
/// sequences unless [`alphabet`](NJConfig::alphabet) is explicitly set;
/// [`substitution_model`](NJConfig::substitution_model) must be compatible with
/// the alphabet or an error is returned.
#[derive(Serialize, Deserialize, ts_rs::TS, Clone, Debug)]
#[ts(export, export_to = "../../wasm/types/lib_types.ts")]
pub struct NJConfig {
    /// The aligned sequences to build the tree from. All sequences must have
    /// the same length.
    pub msa: Vec<SequenceObject>,
    /// Number of bootstrap replicates used to compute support values on
    /// internal nodes. Set to `0` to skip bootstrapping entirely.
    pub n_bootstrap_samples: usize,
    /// Substitution model used to compute pairwise distances.
    pub substitution_model: SubstitutionModel,
    /// Override automatic alphabet detection. When `None` (the default), the
    /// alphabet is inferred from the sequences. Set to `Some(Alphabet::DNA)` or
    /// `Some(Alphabet::Protein)` to force a specific alphabet.
    #[serde(default)]
    pub alphabet: Option<Alphabet>,
    /// Maximum number of threads to use for parallel bootstrap and distance computation.
    ///
    /// Only effective when the `parallel` Cargo feature is enabled.
    /// When `None` (the default), Rayon uses all available hardware threads.
    #[serde(default)]
    pub num_threads: Option<usize>,
}
