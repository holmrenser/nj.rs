//! Multiple sequence alignment container and bootstrap resampling.
//!
//! [`MSA<A>`] stores sequences pre-encoded into typed symbol vectors so that
//! model calculations never need to re-parse raw characters. The type parameter
//! `A` binds the alphabet at compile time, ensuring model–alphabet compatibility
//! is enforced without runtime branching inside the distance computation.

use crate::DistMat;
use crate::alphabet::AlphabetEncoding;
use crate::models::ModelCalculation;
use nanorand::{Rng, WyRand};
use std::collections::HashMap;

/// A typed multiple sequence alignment (MSA).
///
/// Sequences are stored as pre-encoded symbol vectors rather than raw strings,
/// so model calculations never need to re-parse characters. All sequences must
/// have the same length (ragged alignments are rejected at insertion time).
///
/// The type parameter `A` selects the alphabet; use `MSA::<DNA>` or
/// `MSA::<Protein>`. Build an MSA from raw strings with [`push`](MSA::push),
/// [`from_unnamed_sequences`](MSA::from_unnamed_sequences), or by collecting
/// from an iterator of `(String, String)` tuples via the [`FromIterator`] impl.
pub struct MSA<A: AlphabetEncoding> {
    /// Sequence identifiers in insertion order.
    pub identifiers: Vec<String>,
    /// Encoded sequences in insertion order; each inner `Vec` has length
    /// `n_characters`.
    pub sequences: Vec<Vec<A::Symbol>>,
    /// Number of sequences currently in the alignment.
    pub n_sequences: usize,
    /// Alignment length (number of columns), i.e. the length of each sequence.
    pub n_characters: usize,
}

impl<A: AlphabetEncoding> MSA<A> {
    /// Creates a new, empty MSA.
    pub fn new() -> Self {
        Self {
            identifiers: Vec::new(),
            sequences: Vec::new(),
            n_sequences: 0,
            n_characters: 0,
        }
    }

    /// Checks if the MSA is empty.
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    /// Returns the number of sequences in the MSA.
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    /// Returns `Err` if the MSA is non-empty and `seq_len` differs from `n_characters`.
    fn validate_sequence_length(&self, seq_len: usize) -> Result<(), String> {
        if self.n_sequences > 0 && seq_len != self.n_characters {
            return Err(format!(
                "All sequences must have the same length. Expected {}, got {}",
                self.n_characters, seq_len
            ));
        }

        Ok(())
    }

    /// Encodes `seq` and appends it to the MSA. Returns `Err` if `seq` has a
    /// different length than the sequences already present.
    pub fn push(&mut self, id: String, seq: String) -> Result<(), String> {
        let encoded: Vec<A::Symbol> = seq.as_bytes().iter().map(|&b| A::encode(b)).collect();
        self.validate_sequence_length(encoded.len())?;
        self.identifiers.push(id);
        self.sequences.push(encoded);
        self.n_sequences += 1;
        self.n_characters = seq.len();
        Ok(())
    }

    /// Appends an already-encoded sequence. Used by [`bootstrap`](MSA::bootstrap)
    /// to avoid re-encoding symbols that were derived from existing MSA data.
    fn push_encoded(&mut self, id: String, seq: Vec<A::Symbol>) -> Result<(), String> {
        let seq_len = seq.len();
        self.validate_sequence_length(seq_len)?;
        self.identifiers.push(id);
        self.sequences.push(seq);
        self.n_sequences += 1;
        self.n_characters = seq_len;
        Ok(())
    }

    /// Creates an MSA from raw sequences, assigning `Seq0`, `Seq1`, … as identifiers.
    pub fn from_unnamed_sequences(sequences: Vec<String>) -> Result<Self, String> {
        let mut msa = MSA::new();
        for (i, seq) in sequences.into_iter().enumerate() {
            msa.push(format!("Seq{}", i), seq)?;
        }
        Ok(msa)
    }

    /// Returns a map from sequence identifier to its 0-based position in the MSA.
    pub fn to_index_map(&self) -> HashMap<String, usize> {
        self.identifiers
            .iter()
            .enumerate()
            .map(|(i, identifier)| (identifier.clone(), i))
            .collect()
    }

    /// Generates a bootstrap replicate by resampling alignment columns with replacement.
    ///
    /// Draws `n_characters` column indices uniformly at random (with replacement)
    /// using a cryptographically seeded PRNG (`getrandom` → `WyRand`). The
    /// `getrandom` backend is WASM-compatible (uses the JS `crypto` API in
    /// that environment). Returns `Err` if entropy collection fails.
    pub fn bootstrap(&self) -> Result<Self, String> {
        let mut seed_bytes = [0u8; 8];
        getrandom::getrandom(&mut seed_bytes).map_err(|e| format!("RNG error: {e}"))?;
        let mut rng = WyRand::new_seed(u64::from_le_bytes(seed_bytes));
        let sampled_indices: Vec<usize> = (0..self.n_characters)
            .map(|_| rng.generate_range(0..self.n_characters))
            .collect();

        let mut new_msa = MSA::new();
        for (orig_seq, identifier) in self.sequences.iter().zip(self.identifiers.iter()) {
            let new_sequence = sampled_indices
                .iter()
                .map(|&i| *orig_seq.iter().nth(i).unwrap())
                .collect();
            new_msa.push_encoded(identifier.clone(), new_sequence)?;
        }
        Ok(new_msa)
    }

    /// Computes a pairwise distance matrix using substitution model `M`.
    ///
    /// Convenience wrapper around [`DistMat::from_msa`]. Choose `M` to match
    /// the alphabet: DNA models (`JukesCantor`, `Kimura2P`, `PDiff`) for
    /// `MSA::<DNA>`, and protein models (`Poisson`, `PDiff`) for
    /// `MSA::<Protein>`.
    pub fn into_dist<M>(&self) -> DistMat
    where
        M: ModelCalculation<A>,
    {
        DistMat::from_msa::<M, A>(self)
    }
}

impl<A: AlphabetEncoding> IntoIterator for MSA<A> {
    type Item = (String, Vec<A::Symbol>);
    type IntoIter = std::iter::Zip<std::vec::IntoIter<String>, std::vec::IntoIter<Vec<A::Symbol>>>;

    fn into_iter(self) -> Self::IntoIter {
        self.identifiers.into_iter().zip(self.sequences.into_iter())
    }
}

impl<A: AlphabetEncoding> FromIterator<(String, String)> for MSA<A> {
    fn from_iter<I: IntoIterator<Item = (String, String)>>(iter: I) -> Self {
        let mut msa = MSA::new();
        for (id, seq) in iter {
            msa.push(id, seq)
                .expect("All sequences in an MSA must have the same length");
        }
        msa
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::{DNA, Protein};

    #[test]
    fn test_msa_push_and_len() {
        let mut msa = MSA::<DNA>::new();
        msa.push("seq1".into(), "ACGT".into()).unwrap();
        msa.push("seq2".into(), "AGGT".into()).unwrap();
        assert_eq!(msa.len(), 2);
        assert_eq!(msa.n_characters, 4);
    }

    #[test]
    fn test_msa_push_rejects_ragged_alignment() {
        let mut msa = MSA::<DNA>::new();
        msa.push("seq1".into(), "ACGT".into()).unwrap();
        let err = msa.push("seq2".into(), "AGG".into()).unwrap_err();
        assert!(err.contains("All sequences must have the same length"));
    }

    #[test]
    fn test_msa_bootstrap() {
        let mut msa = MSA::<DNA>::new();
        msa.push("seq1".into(), "ACGT".into()).unwrap();
        msa.push("seq2".into(), "AGGT".into()).unwrap();

        let boot_msa = msa.bootstrap().unwrap();
        assert_eq!(boot_msa.len(), 2);
        assert_eq!(boot_msa.n_characters, 4);
    }

    #[test]
    fn test_msa_to_index_map() {
        let mut msa = MSA::<DNA>::new();
        msa.push("seq1".into(), "ACGT".into()).unwrap();
        msa.push("seq2".into(), "AGGT".into()).unwrap();
        let index_map = msa.to_index_map();
        assert_eq!(index_map.get("seq1"), Some(&0));
        assert_eq!(index_map.get("seq2"), Some(&1));
    }

    #[test]
    fn test_msa_from_unnamed_sequences() {
        let sequences = vec!["ACGT".into(), "AGGT".into()];
        let msa = MSA::<DNA>::from_unnamed_sequences(sequences).unwrap();
        assert_eq!(msa.len(), 2);
        assert_eq!(msa.identifiers[0], "Seq0");
        assert_eq!(msa.identifiers[1], "Seq1");
    }

    #[test]
    fn test_msa_into_iterator() {
        let mut msa = MSA::<DNA>::new();
        msa.push("seq1".into(), "ACGT".into()).unwrap();
        msa.push("seq2".into(), "AGGT".into()).unwrap();
        let mut iter = msa.into_iter();
        let (id1, _) = iter.next().unwrap();
        assert_eq!(id1, "seq1");
        let (id2, _) = iter.next().unwrap();
        assert_eq!(id2, "seq2");
    }

    #[test]
    fn test_msa_from_iterator() {
        let data = vec![
            ("seq1".into(), "ACGT".into()),
            ("seq2".into(), "AGGT".into()),
        ];
        let msa: MSA<DNA> = data.into_iter().collect();
        assert_eq!(msa.len(), 2);
        assert_eq!(msa.identifiers[0], "seq1");
        assert_eq!(msa.identifiers[1], "seq2");
    }

    #[test]
    fn test_msa_is_empty() {
        let msa = MSA::<DNA>::new();
        assert!(msa.is_empty());
        let mut msa2 = MSA::<DNA>::new();
        msa2.push("seq1".into(), "ACGT".into()).unwrap();
        assert!(!msa2.is_empty());
    }

    #[test]
    fn test_msa_protein() {
        let mut msa = MSA::<Protein>::new();
        msa.push("prot1".into(), "ACDEFGHIKLMNPQRSTVWY".into())
            .unwrap();
        msa.push("prot2".into(), "ACDEFGHIKLMNPQRSTVWY".into())
            .unwrap();
        assert_eq!(msa.len(), 2);
        assert_eq!(msa.n_characters, 20);
    }

    #[test]
    fn test_msa_bootstrap_protein() {
        let mut msa = MSA::<Protein>::new();
        msa.push("prot1".into(), "ACDEFGHIKLMNPQRSTVWY".into())
            .unwrap();
        msa.push("prot2".into(), "ACDEFGHIKLMNPQRSTVWY".into())
            .unwrap();
        let boot_msa = msa.bootstrap().unwrap();
        assert_eq!(boot_msa.len(), 2);
        assert_eq!(boot_msa.n_characters, 20);
    }

    #[test]
    fn test_msa_to_index_map_protein() {
        let mut msa = MSA::<Protein>::new();
        msa.push("prot1".into(), "ACDEFGHIKLMNPQRSTVWY".into())
            .unwrap();
        msa.push("prot2".into(), "ACDEFGHIKLMNPQRSTVWY".into())
            .unwrap();
        let index_map = msa.to_index_map();
        assert_eq!(index_map.get("prot1"), Some(&0));
        assert_eq!(index_map.get("prot2"), Some(&1));
    }

    #[test]
    fn test_msa_from_unnamed_sequences_protein() {
        let sequences = vec!["ACDEFGHIKLMNPQRSTVWY".into(), "ACDEFGHIKLMNPQRSTVWY".into()];
        let msa = MSA::<Protein>::from_unnamed_sequences(sequences).unwrap();
        assert_eq!(msa.len(), 2);
        assert_eq!(msa.identifiers[0], "Seq0");
        assert_eq!(msa.identifiers[1], "Seq1");
    }

    #[test]
    fn test_msa_into_iterator_protein() {
        let mut msa = MSA::<Protein>::new();
        msa.push("prot1".into(), "ACDEFGHIKLMNPQRSTVWY".into())
            .unwrap();
        msa.push("prot2".into(), "ACDEFGHIKLMNPQRSTVWY".into())
            .unwrap();
        let mut iter = msa.into_iter();
        let (id1, _) = iter.next().unwrap();
        assert_eq!(id1, "prot1");
        let (id2, _) = iter.next().unwrap();
        assert_eq!(id2, "prot2");
    }
    #[test]
    fn test_msa_from_iterator_protein() {
        let data = vec![
            ("prot1".into(), "ACDEFGHIKLMNPQRSTVWY".into()),
            ("prot2".into(), "ACDEFGHIKLMNPQRSTVWY".into()),
        ];
        let msa: MSA<Protein> = data.into_iter().collect();
        assert_eq!(msa.len(), 2);
        assert_eq!(msa.identifiers[0], "prot1");
        assert_eq!(msa.identifiers[1], "prot2");
    }

    #[test]
    fn test_msa_is_empty_protein() {
        let msa = MSA::<Protein>::new();
        assert!(msa.is_empty());
        let mut msa2 = MSA::<Protein>::new();
        msa2.push("prot1".into(), "ACDEFGHIKLMNPQRSTVWY".into())
            .unwrap();
        assert!(!msa2.is_empty());
    }

    #[test]
    fn test_msa_into_dist_dna() {
        let mut msa = MSA::<DNA>::new();
        msa.push("seq1".into(), "ACGT".into()).unwrap();
        msa.push("seq2".into(), "AGGT".into()).unwrap();
        let dist_mat = msa.into_dist::<crate::models::PDiff>();
        assert_eq!(dist_mat.names.len(), 2);
    }

    #[test]
    fn test_msa_into_dist_protein() {
        let mut msa = MSA::<Protein>::new();
        msa.push("prot1".into(), "ACDEFGHIKLMNPQRSTVWY".into())
            .unwrap();
        msa.push("prot2".into(), "ACDEFGHIKLMNPQRSTVWY".into())
            .unwrap();
        let dist_mat = msa.into_dist::<crate::models::PDiff>();
        assert_eq!(dist_mat.names.len(), 2);
    }
}
