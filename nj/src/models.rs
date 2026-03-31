use serde::{Deserialize, Serialize};

use crate::alphabet::{AlphabetEncoding, DNA, DnaSymbol, Protein, ProteinSymbol};

pub trait ModelCalculation<A: AlphabetEncoding> {
    type Acc;

    fn init() -> Self::Acc;
    fn accumulate(acc: &mut Self::Acc, a: A::Symbol, b: A::Symbol) -> Self::Acc;
    fn finalize(acc: &Self::Acc, aln_len: usize) -> f64;
}

#[inline(always)]
pub fn pairwise_distance<M, A>(s1: &[A::Symbol], s2: &[A::Symbol]) -> f64
where
    M: ModelCalculation<A>,
    A: AlphabetEncoding,
{
    let aln_len = s1.len();
    let mut acc = M::init();

    for k in 0..aln_len {
        acc = M::accumulate(&mut acc, s1[k], s2[k]);
    }

    M::finalize(&acc, aln_len)
}

pub struct PDiff;

impl ModelCalculation<DNA> for PDiff {
    type Acc = usize; // count of differences

    fn init() -> Self::Acc {
        0
    }

    fn accumulate(acc: &mut Self::Acc, a: DnaSymbol, b: DnaSymbol) -> Self::Acc {
        if a != b && a != DnaSymbol::Gap && b != DnaSymbol::Gap {
            *acc += 1;
        }
        *acc
    }

    fn finalize(acc: &Self::Acc, aln_len: usize) -> f64 {
        *acc as f64 / aln_len as f64
    }
}

impl ModelCalculation<Protein> for PDiff {
    type Acc = usize; // count of differences

    fn init() -> Self::Acc {
        0
    }

    fn accumulate(acc: &mut Self::Acc, a: ProteinSymbol, b: ProteinSymbol) -> Self::Acc {
        if a != b && a != ProteinSymbol::Gap && b != ProteinSymbol::Gap {
            *acc += 1;
        }
        *acc
    }

    fn finalize(acc: &Self::Acc, aln_len: usize) -> f64 {
        *acc as f64 / aln_len as f64
    }
}

pub struct JukesCantor;

impl ModelCalculation<DNA> for JukesCantor {
    type Acc = usize; // count of differences

    fn init() -> Self::Acc {
        0
    }

    fn accumulate(acc: &mut Self::Acc, a: DnaSymbol, b: DnaSymbol) -> Self::Acc {
        if a != b && a != DnaSymbol::Gap && b != DnaSymbol::Gap {
            *acc += 1;
        }
        *acc
    }

    fn finalize(acc: &Self::Acc, aln_len: usize) -> f64 {
        let p = *acc as f64 / aln_len as f64;
        if p >= 0.75 {
            f64::INFINITY // distance undefined
        } else {
            -0.75 * (1.0 - (4.0 / 3.0) * p).ln()
        }
    }
}

pub struct Kimura2P;

impl ModelCalculation<DNA> for Kimura2P {
    type Acc = (usize, usize); // (transitions, transversions)

    fn init() -> Self::Acc {
        (0, 0)
    }

    fn accumulate(acc: &mut Self::Acc, a: DnaSymbol, b: DnaSymbol) -> Self::Acc {
        if a != b && a != DnaSymbol::Gap && b != DnaSymbol::Gap {
            match (a, b) {
                (DnaSymbol::A, DnaSymbol::G)
                | (DnaSymbol::G, DnaSymbol::A)
                | (DnaSymbol::C, DnaSymbol::T)
                | (DnaSymbol::T, DnaSymbol::C) => acc.0 += 1, // transition
                _ => acc.1 += 1, // transversion
            }
        }
        *acc
    }

    fn finalize(acc: &Self::Acc, aln_len: usize) -> f64 {
        let (ti, tv) = *acc;
        let p = ti as f64 / aln_len as f64;
        let q = tv as f64 / aln_len as f64;
        let denom1 = 1.0 - 2.0 * p - q;
        let denom2 = 1.0 - 2.0 * q;
        if denom1 <= 0.0 || denom2 <= 0.0 {
            f64::INFINITY // distance undefined
        } else {
            -0.5 * denom1.ln() - 0.25 * denom2.ln()
        }
    }
}

pub struct Poisson;

impl ModelCalculation<Protein> for Poisson {
    type Acc = usize; // count of differences

    fn init() -> Self::Acc {
        0
    }

    fn accumulate(acc: &mut Self::Acc, a: ProteinSymbol, b: ProteinSymbol) -> Self::Acc {
        if a != b && a != ProteinSymbol::Gap && b != ProteinSymbol::Gap {
            *acc += 1;
        }
        *acc
    }

    fn finalize(acc: &Self::Acc, aln_len: usize) -> f64 {
        let p = *acc as f64 / aln_len as f64;
        if p >= 1.0 {
            f64::INFINITY // distance undefined
        } else {
            -(1.0 - p).ln()
        }
    }
}

#[derive(Clone, Debug, ts_rs::TS, Serialize, Deserialize)]
#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
pub enum SubstitutionModel {
    PDiff,
    JukesCantor,
    Kimura2P,
    Poisson,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::{DNA, Protein};

    #[test]
    fn test_pdiff_dna() {
        let s1 = vec![
            DnaSymbol::A,
            DnaSymbol::C,
            DnaSymbol::G,
            DnaSymbol::T,
            DnaSymbol::A,
        ];
        let s2 = vec![
            DnaSymbol::A,
            DnaSymbol::G,
            DnaSymbol::G,
            DnaSymbol::T,
            DnaSymbol::C,
        ];
        let dist = pairwise_distance::<PDiff, DNA>(&s1, &s2);
        assert_eq!(dist, 0.4); // 2 differences out of 5
    }

    #[test]
    fn test_jukes_cantor_dna() {
        let s1 = vec![
            DnaSymbol::A,
            DnaSymbol::C,
            DnaSymbol::G,
            DnaSymbol::T,
            DnaSymbol::A,
        ];
        let s2 = vec![
            DnaSymbol::A,
            DnaSymbol::G,
            DnaSymbol::G,
            DnaSymbol::T,
            DnaSymbol::C,
        ];
        let dist = pairwise_distance::<JukesCantor, DNA>(&s1, &s2);
        let p = 0.4;
        let expected = -0.75 * ((1.0 - (4.0 / 3.0) * p) as f64).ln();
        assert!((dist - expected).abs() < 1e-6);
    }
}
