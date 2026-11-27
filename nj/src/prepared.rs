//! A prepared alignment: the encoded sequences plus the per-sequence data that makes
//! pair counting fast.
//!
//! For DNA, every sequence is also stored as four **bit-planes** (A, C, G, T): bit `i`
//! of plane `a` is set iff column `i` holds base `a`. Columns with a gap, `N` or an
//! ambiguity code set no bit, so they drop out of every count automatically (pairwise
//! deletion). The pair-count table is then
//!
//! ```text
//! F[a][b] = Σ_words popcount(plane_a(i) & plane_b(j))
//! ```
//!
//! which handles 64 columns per AND and popcount, instead of one read-modify-write on
//! a data-dependent histogram cell per column.

use crate::alphabet::{AlphabetEncoding, DNA};
use crate::counts::DnaCounts;
use crate::kernels;
use crate::msa::MSA;
use fearless_simd::Level;

/// Bit-planes of every sequence of a DNA alignment.
struct Planes {
    /// 64-bit words per plane.
    words: usize,
    /// `data[(seq * 4 + base) * words + k]`: word `k` of plane `base` of sequence `seq`.
    data: Vec<u64>,
}

impl Planes {
    fn build<'s>(sequences: impl Iterator<Item = &'s [u8]>, n_sequences: usize, n_characters: usize) -> Self {
        let words = n_characters.div_ceil(64);
        let mut data = vec![0u64; n_sequences * 4 * words];
        for (s, bytes) in sequences.enumerate() {
            let planes = &mut data[s * 4 * words..(s + 1) * 4 * words];
            for (i, &b) in bytes.iter().enumerate() {
                if b < 4 {
                    planes[b as usize * words + i / 64] |= 1u64 << (i % 64);
                }
            }
        }
        Planes { words, data }
    }

    /// The four planes of sequence `s`, concatenated.
    #[inline]
    fn sequence(&self, s: usize) -> &[u64] {
        &self.data[s * 4 * self.words..(s + 1) * 4 * self.words]
    }
}

/// An alignment prepared for distance computation. Borrows the encoded sequences of
/// an [`MSA`]; for DNA it additionally owns the bit-planes.
pub struct Prepared<'a, A: AlphabetEncoding> {
    sequences: &'a [Vec<A::Symbol>],
    planes: Option<Planes>,
    /// SIMD level chosen once, so pair counting does not detect it per pair.
    level: Level,
}

impl<'a, A: AlphabetEncoding> Prepared<'a, A> {
    /// Prepares `msa`. Builds bit-planes when the alphabet has four states (DNA).
    pub fn new(msa: &'a MSA<A>) -> Self {
        Self::new_at(msa, Level::new())
    }

    /// Like [`new`](Prepared::new) with an explicit SIMD level (tests run every level).
    pub(crate) fn new_at(msa: &'a MSA<A>, level: Level) -> Self {
        let planes = (A::N_STATES == 4).then(|| {
            Planes::build(
                msa.sequences.iter().map(|s| A::as_bytes(s)),
                msa.sequences.len(),
                msa.n_characters,
            )
        });
        Prepared {
            sequences: &msa.sequences,
            planes,
            level,
        }
    }

    /// Number of sequences.
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    /// `true` if the alignment has no sequences.
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    /// The encoded symbols of sequence `i`.
    #[inline]
    pub fn sequence(&self, i: usize) -> &'a [A::Symbol] {
        &self.sequences[i]
    }
}

impl Prepared<'_, DNA> {
    /// The pair-count table of sequences `i` and `j`, from the bit-planes.
    #[inline]
    pub fn dna_counts(&self, i: usize, j: usize) -> DnaCounts {
        let planes = self.planes.as_ref().expect("DNA alignments always have planes");
        DnaCounts {
            f: kernels::plane_counts(self.level, planes.sequence(i), planes.sequence(j), planes.words),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
    }

    fn random_dna(rng: &mut Rng, n: usize, len: usize) -> MSA<DNA> {
        let pairs = (0..n).map(|i| {
            let seq: String = (0..len)
                .map(|_| {
                    // mostly bases, some N, gaps and IUPAC codes
                    match rng.next() % 20 {
                        0 => 'N',
                        1 => '-',
                        2 => 'R',
                        r => b"ACGT"[(r % 4) as usize] as char,
                    }
                })
                .collect();
            (format!("s{i}"), seq)
        });
        MSA::<DNA>::try_from_pairs(pairs).unwrap()
    }

    #[test]
    fn test_plane_counts_equal_the_histogram() {
        let mut rng = Rng(0x1234_5678_9ABC_DEF1);
        // Lengths around the word boundaries, and a few long ones.
        for len in [1usize, 2, 7, 63, 64, 65, 127, 128, 129, 191, 193, 640, 1000] {
            let msa = random_dna(&mut rng, 6, len);
            for (name, level) in crate::kernels::testing::levels() {
                let prepared = Prepared::new_at(&msa, level);
                for i in 0..6 {
                    for j in 0..6 {
                        let want = DnaCounts::from_symbols(&msa.sequences[i], &msa.sequences[j]);
                        assert_eq!(prepared.dna_counts(i, j), want, "level {name}, len {len}, pair ({i},{j})");
                    }
                }
            }
        }
    }

    #[test]
    fn test_prepared_protein_has_no_planes() {
        let msa = MSA::<crate::alphabet::Protein>::try_from_pairs(
            [("a".to_string(), "ARND".to_string()), ("b".to_string(), "ARNE".to_string())].into_iter(),
        )
        .unwrap();
        let prepared = Prepared::new(&msa);
        assert!(prepared.planes.is_none());
        assert_eq!(prepared.len(), 2);
        assert_eq!(prepared.sequence(1).len(), 4);
    }
}
