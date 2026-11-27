//! Input validation and alphabet detection.

use std::collections::HashSet;

use crate::alphabet::Alphabet;
use crate::config::{DEFAULT_MAX_DISTANCE, SequenceObject};
use crate::error::NJError;
use crate::event::{LogLevel, NJEvent};
use crate::models::{RateHet, SubstitutionModel};

/// The per-model options of a run: the rate matrix of `Custom`, the exchangeabilities and base
/// frequencies of `Gtr`. Borrowed from the configuration.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct ModelOptions<'a> {
    pub matrix: Option<&'a str>,
    pub gtr_rates: Option<&'a [f64]>,
    pub base_frequencies: Option<&'a [f64]>,
}

/// Checks that each option is present exactly for the model that uses it: `matrix` for `Custom`,
/// `gtr_rates` for `Gtr` (`base_frequencies` is optional there), and none for any other model.
pub(crate) fn validate_model_options(model: &SubstitutionModel, opts: &ModelOptions<'_>) -> Result<(), NJError> {
    let err = |reason: &str| Err(NJError::InvalidMatrix { reason: reason.to_string() });
    match model {
        SubstitutionModel::Custom => {
            if opts.matrix.is_none() {
                return err("the 'Custom' model needs a rate matrix in the `matrix` option");
            }
        }
        SubstitutionModel::Gtr => {
            if opts.gtr_rates.is_none() {
                return err("the 'Gtr' model needs its six exchangeabilities (AC, AG, AT, CG, CT, GT) in `gtr_rates`");
            }
        }
        _ => {}
    }
    if opts.matrix.is_some() && !matches!(model, SubstitutionModel::Custom) {
        return err("`matrix` is only used with the 'Custom' model");
    }
    if opts.gtr_rates.is_some() && !matches!(model, SubstitutionModel::Gtr) {
        return err("`gtr_rates` is only used with the 'Gtr' model");
    }
    if opts.base_frequencies.is_some() && !matches!(model, SubstitutionModel::Gtr) {
        return err("`base_frequencies` is only used with the 'Gtr' model");
    }
    Ok(())
}

/// Validates that the MSA is non-empty, all sequences have equal length, and
/// all identifiers are unique.
pub(crate) fn validate_msa(msa: &[SequenceObject]) -> Result<(), NJError> {
    if msa.is_empty() {
        return Err(NJError::EmptyMsa);
    }
    let expected_len = msa[0].sequence.len();
    if expected_len == 0 {
        return Err(NJError::EmptySequence);
    }
    let mut seen: HashSet<&str> = HashSet::with_capacity(msa.len());
    for s in msa {
        if s.sequence.len() != expected_len {
            return Err(NJError::SequenceLengthMismatch {
                expected: expected_len,
                got: s.sequence.len(),
                identifier: s.identifier.clone(),
            });
        }
        if !seen.insert(s.identifier.as_str()) {
            return Err(NJError::DuplicateIdentifier {
                identifier: s.identifier.clone(),
            });
        }
    }
    Ok(())
}

/// Validates the among-site rate-variation parameters and packs them into a
/// [`RateHet`].
///
/// `gamma_shape`, if present, must be finite and strictly positive; `p_invar`,
/// if present, must be finite and in `[0, 1)`. A missing `p_invar` defaults to
/// `0.0` (no invariant-sites correction).
pub(crate) fn validate_rate_params(
    gamma_shape: Option<f64>,
    p_invar: Option<f64>,
) -> Result<RateHet, NJError> {
    if let Some(alpha) = gamma_shape
        && (!alpha.is_finite() || alpha <= 0.0)
    {
        return Err(NJError::InvalidGammaShape { value: alpha });
    }
    if let Some(pinv) = p_invar
        && (!pinv.is_finite() || !(0.0..1.0).contains(&pinv))
    {
        return Err(NJError::InvalidPInvar { value: pinv });
    }
    Ok(RateHet {
        gamma_shape,
        p_invar: p_invar.unwrap_or(0.0),
    })
}

/// Validates `max_distance` (finite and positive) and applies the default when absent.
pub(crate) fn validate_max_distance(max_distance: Option<f64>) -> Result<f64, NJError> {
    match max_distance {
        None => Ok(DEFAULT_MAX_DISTANCE),
        Some(d) if d.is_finite() && d > 0.0 => Ok(d),
        Some(d) => Err(NJError::InvalidMaxDistance { value: d }),
    }
}

/// Heuristically detects whether the MSA contains DNA or protein sequences.
///
/// Returns [`Alphabet::DNA`] unless any sequence contains a byte that is not
/// in the DNA character set (case-insensitive), in which case
/// [`Alphabet::Protein`] is returned. The DNA set includes standard bases
/// `{A, C, G, T}`, uridine `U` (RNA), `N` (unknown), gap `-`, and all 11
/// IUPAC ambiguity codes `{R, Y, S, W, K, M, B, D, H, V}`.
pub(crate) fn detect_alphabet(msa: &[SequenceObject]) -> Alphabet {
    let mut is_protein = false;

    'outer: for seq in msa {
        for c in seq.sequence.bytes() {
            match c.to_ascii_uppercase() {
                b'A' | b'C' | b'G' | b'T' | b'U' | b'N' | b'-' | b'R' | b'Y' | b'S' | b'W'
                | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' => { /* still possible DNA */ }
                _ => {
                    is_protein = true;
                    break 'outer;
                }
            }
        }
    }

    if is_protein {
        Alphabet::Protein
    } else {
        Alphabet::DNA
    }
}

/// Emits a `Warning`-level [`NJEvent::Log`] for any sequence consisting entirely
/// of gap characters. Such sequences have no comparable sites, so every pairwise
/// distance involving them is `0.0` ([`crate::models`]); surfacing a warning makes
/// that otherwise-silent behaviour observable. No-op when `on_event` is `None`.
pub(crate) fn warn_degenerate_sequences(
    msa: &[SequenceObject],
    on_event: Option<&dyn Fn(NJEvent)>,
) {
    let Some(cb) = on_event else { return };
    for s in msa {
        if !s.sequence.is_empty() && s.sequence.bytes().all(|b| b == b'-') {
            cb(NJEvent::Log {
                level: LogLevel::Warning,
                message: format!(
                    "Sequence '{}' contains only gaps; its pairwise distances will be 0.0",
                    s.identifier
                ),
            });
        }
    }
}
