use crate::DistMat;
use crate::alphabet::AlphabetEncoding;
use crate::models::ModelCalculation;
use nanorand::Rng;
use nanorand::WyRand;
use std::collections::HashMap;

pub struct MSA<A: AlphabetEncoding> {
    pub identifiers: Vec<String>,
    pub sequences: Vec<Vec<A::Symbol>>,
    pub n_sequences: usize,
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

    /// Adds a FastaSequence to the MSA.
    pub fn push(&mut self, id: String, seq: String) {
        let encoded: Vec<A::Symbol> = seq.as_bytes().iter().map(|&b| A::encode(b)).collect();
        self.identifiers.push(id);
        self.sequences.push(encoded);
        self.n_sequences += 1;
        self.n_characters = seq.len();
    }

    fn push_encoded(&mut self, id: String, seq: Vec<A::Symbol>) {
        self.identifiers.push(id);
        self.sequences.push(seq);
        self.n_sequences += 1;
        self.n_characters = self.sequences[0].len();
    }

    /// Creates an MSA from a vector of unnamed sequences, assigning default names.
    pub fn from_unnamed_sequences(sequences: Vec<String>) -> Self {
        let mut msa = MSA::new();
        for (i, seq) in sequences.into_iter().enumerate() {
            msa.push(format!("Seq{}", i), seq);
        }
        msa
    }

    /// Creates an index map from sequence identifiers to their indices in the MSA.
    pub fn to_index_map(&self) -> HashMap<String, usize> {
        self.identifiers
            .iter()
            .enumerate()
            .map(|(i, identifier)| (identifier.clone(), i))
            .collect()
    }

    /// Generates a bootstrap replicate of the MSA by resampling columns with replacement.
    pub fn bootstrap(&self) -> Self {
        let mut rng = WyRand::new();
        let sampled_indices: Vec<usize> = (0..self.n_characters)
            .map(|_| rng.generate_range(0..self.n_characters))
            .collect();

        let mut new_msa = MSA::new();
        for (orig_seq, identifier) in self.sequences.iter().zip(self.identifiers.iter()) {
            let new_sequence = sampled_indices
                .iter()
                .map(|&i| *orig_seq.iter().nth(i).unwrap())
                .collect();
            new_msa.push_encoded(identifier.clone(), new_sequence);
        }
        new_msa
    }

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
            msa.push(id, seq);
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
        msa.push("seq1".into(), "ACGT".into());
        msa.push("seq2".into(), "AGGT".into());
        assert_eq!(msa.len(), 2);
        assert_eq!(msa.n_characters, 4);
    }

    #[test]
    fn test_msa_bootstrap() {
        let mut msa = MSA::<DNA>::new();
        msa.push("seq1".into(), "ACGT".into());
        msa.push("seq2".into(), "AGGT".into());

        let boot_msa = msa.bootstrap();
        assert_eq!(boot_msa.len(), 2);
        assert_eq!(boot_msa.n_characters, 4);
    }

    #[test]
    fn test_msa_to_index_map() {
        let mut msa = MSA::<DNA>::new();
        msa.push("seq1".into(), "ACGT".into());
        msa.push("seq2".into(), "AGGT".into());
        let index_map = msa.to_index_map();
        assert_eq!(index_map.get("seq1"), Some(&0));
        assert_eq!(index_map.get("seq2"), Some(&1));
    }

    #[test]
    fn test_msa_from_unnamed_sequences() {
        let sequences = vec!["ACGT".into(), "AGGT".into()];
        let msa = MSA::<DNA>::from_unnamed_sequences(sequences);
        assert_eq!(msa.len(), 2);
        assert_eq!(msa.identifiers[0], "Seq0");
        assert_eq!(msa.identifiers[1], "Seq1");
    }

    #[test]
    fn test_msa_into_iterator() {
        let mut msa = MSA::<DNA>::new();
        msa.push("seq1".into(), "ACGT".into());
        msa.push("seq2".into(), "AGGT".into());
        let mut iter = msa.into_iter();
        let (id1, seq1) = iter.next().unwrap();
        assert_eq!(id1, "seq1");
        let (id2, seq2) = iter.next().unwrap();
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
        msa2.push("seq1".into(), "ACGT".into());
        assert!(!msa2.is_empty());
    }

    #[test]
    fn test_msa_protein() {
        let mut msa = MSA::<Protein>::new();
        msa.push("prot1".into(), "ACDEFGHIKLMNPQRSTVWY".into());
        msa.push("prot2".into(), "ACDEFGHIKLMNPQRSTVWY".into());
        assert_eq!(msa.len(), 2);
        assert_eq!(msa.n_characters, 20);
    }

    #[test]
    fn test_msa_bootstrap_protein() {
        let mut msa = MSA::<Protein>::new();
        msa.push("prot1".into(), "ACDEFGHIKLMNPQRSTVWY".into());
        msa.push("prot2".into(), "ACDEFGHIKLMNPQRSTVWY".into());
        let boot_msa = msa.bootstrap();
        assert_eq!(boot_msa.len(), 2);
        assert_eq!(boot_msa.n_characters, 20);
    }

    #[test]
    fn test_msa_to_index_map_protein() {
        let mut msa = MSA::<Protein>::new();
        msa.push("prot1".into(), "ACDEFGHIKLMNPQRSTVWY".into());
        msa.push("prot2".into(), "ACDEFGHIKLMNPQRSTVWY".into());
        let index_map = msa.to_index_map();
        assert_eq!(index_map.get("prot1"), Some(&0));
        assert_eq!(index_map.get("prot2"), Some(&1));
    }

    #[test]
    fn test_msa_from_unnamed_sequences_protein() {
        let sequences = vec!["ACDEFGHIKLMNPQRSTVWY".into(), "ACDEFGHIKLMNPQRSTVWY".into()];
        let msa = MSA::<Protein>::from_unnamed_sequences(sequences);
        assert_eq!(msa.len(), 2);
        assert_eq!(msa.identifiers[0], "Seq0");
        assert_eq!(msa.identifiers[1], "Seq1");
    }

    #[test]
    fn test_msa_into_iterator_protein() {
        let mut msa = MSA::<Protein>::new();
        msa.push("prot1".into(), "ACDEFGHIKLMNPQRSTVW   Y".into());
        msa.push("prot2".into(), "ACDEFGHIKLMNPQRSTVWY".into());
        let mut iter = msa.into_iter();
        let (id1, seq1) = iter.next().unwrap();
        assert_eq!(id1, "prot1");
        let (id2, seq2) = iter.next().unwrap();
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
        msa2.push("prot1".into(), "ACDEFGHIKLMNPQRSTVWY".into());
        assert!(!msa2.is_empty());
    }

    #[test]
    fn test_msa_into_dist_dna() {
        let mut msa = MSA::<DNA>::new();
        msa.push("seq1".into(), "ACGT".into());
        msa.push("seq2".into(), "AGGT".into());
        let dist_mat = msa.into_dist::<crate::models::PDiff>();
        assert_eq!(dist_mat.names.len(), 2);
    }

    #[test]
    fn test_msa_into_dist_protein() {
        let mut msa = MSA::<Protein>::new();
        msa.push("prot1".into(), "ACDEFGHIKLMNPQRSTVWY".into());
        msa.push("prot2".into(), "ACDEFGHIKLMNPQRSTVWY".into());
        let dist_mat = msa.into_dist::<crate::models::PDiff>();
        assert_eq!(dist_mat.names.len(), 2);
    }
}
