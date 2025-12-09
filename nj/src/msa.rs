use crate::DistMat;
use crate::config::FastaSequence;
use nanorand::Rng;
use nanorand::WyRand;

pub struct MSA {
    pub sequences: Vec<FastaSequence>,
    pub n_sequences: usize,
    pub n_characters: usize,
}

impl MSA {
    /// Creates a new, empty MSA.
    pub fn new() -> Self {
        Self {
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
    pub fn push(&mut self, seq: FastaSequence) {
        self.n_sequences += 1;
        self.n_characters = seq.sequence.len();
        self.sequences.push(seq);
    }

    pub fn from_unnamed_sequences(sequences: Vec<String>) -> Self {
        let n_sequences = sequences.len();
        let n_characters = if n_sequences == 0 {
            0
        } else {
            sequences[0].len()
        };
        let s = sequences
            .into_iter()
            .enumerate()
            .map(|(i, s)| FastaSequence {
                identifier: format!("Seq{}", i),
                sequence: s,
            })
            .collect();
        MSA {
            sequences: s,
            n_sequences,
            n_characters,
        }
    }

    /// Generates a bootstrap replicate of the MSA by resampling columns with replacement.
    pub fn bootstrap(&self) -> Self {
        let mut rng = WyRand::new();
        let sampled_indices: Vec<usize> = (0..self.n_characters)
            .map(|_| rng.generate_range(0..self.n_characters))
            .collect();

        let mut new_msa = MSA::new();
        for fs in &self.sequences {
            let new_sequence: String = sampled_indices
                .iter()
                .map(|&i| fs.sequence.chars().nth(i).unwrap())
                .collect();
            new_msa.push(FastaSequence {
                identifier: fs.identifier.clone(),
                sequence: new_sequence,
            });
        }
        new_msa
    }

    /// Converts the MSA to a distance matrix using p-distance (proportion of differing sites).
    pub fn into_dist(&self) -> DistMat {
        let n = self.len();
        let names: Vec<String> = self
            .sequences
            .clone()
            .into_iter()
            .map(|fs| fs.identifier.clone())
            .collect();
        let mut dist = DistMat::empty_with_names(names);

        (0..n).for_each(|i| {
            (0..i).for_each(|j| {
                let (diffs, valid) = self.sequences[i]
                    .sequence
                    .chars()
                    .zip(self.sequences[j].sequence.chars())
                    .filter(|(a, b)| *a != '-' && *b != '-')
                    .fold((0, 0), |(d, v), (a, b)| {
                        (d + if a != b { 1 } else { 0 }, v + 1)
                    });
                let d = if valid > 0 {
                    diffs as f64 / valid as f64
                } else {
                    0.0
                };
                dist.set(i, j, d);
            })
        });
        dist
    }
}

/// Allows iteration over the sequences in the MSA.
impl IntoIterator for MSA {
    type Item = FastaSequence;
    type IntoIter = std::vec::IntoIter<FastaSequence>;

    fn into_iter(self) -> Self::IntoIter {
        self.sequences.into_iter()
    }
}

/// Allows creating an MSA from an iterator of FastaSequence.
impl FromIterator<FastaSequence> for MSA {
    fn from_iter<I: IntoIterator<Item = FastaSequence>>(iter: I) -> Self {
        let mut msa = MSA::new();
        for seq in iter {
            msa.push(seq);
        }
        msa
    }
}

/// Allows indexing into the MSA to get a specific FastaSequence.
impl std::ops::Index<usize> for MSA {
    type Output = FastaSequence;

    fn index(&self, index: usize) -> &Self::Output {
        &self.sequences[index]
    }
}

/// Allows mutable indexing into the MSA to get a specific FastaSequence.
impl std::ops::IndexMut<usize> for MSA {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.sequences[index]
    }
}
