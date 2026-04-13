//! FASTA parser for multiple sequence alignment input.
//!
//! [`parse_fasta`] converts a FASTA-formatted string into a [`Vec<SequenceObject>`]
//! ready to be passed directly to [`crate::nj`], [`crate::distance_matrix`], or
//! [`crate::average_distance`]. The parser is intentionally simple — it handles
//! multi-line sequences, trims whitespace, and skips blank lines, but does not
//! validate sequence characters (that happens inside the library).

use crate::config::SequenceObject;
use crate::error::NJError;

/// Parses a FASTA-formatted string into a vector of [`SequenceObject`]s.
///
/// - Leading and trailing whitespace on each line is stripped.
/// - Multi-line sequences are concatenated.
/// - Blank lines are ignored.
/// - The header text after `>` is used as the sequence identifier (trimmed).
///
/// Returns [`NJError::ParseError`] if:
/// - A sequence line appears before any header.
/// - A header has no sequence (empty or only whitespace after it).
/// - The input contains no sequences at all.
/// - Sequences have different lengths (alignment check).
pub fn parse_fasta(input: &str) -> Result<Vec<SequenceObject>, NJError> {
    let mut msa = Vec::<SequenceObject>::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = String::new();

    for line in input.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(header) = trimmed.strip_prefix('>') {
            // Flush previous entry.
            if let Some(identifier) = current_name.replace(header.trim().to_string()) {
                if current_seq.is_empty() {
                    return Err(NJError::ParseError(format!(
                        "sequence '{identifier}' has no sequence data"
                    )));
                }
                msa.push(SequenceObject { identifier, sequence: current_seq });
                current_seq = String::new();
            }
        } else {
            if current_name.is_none() {
                return Err(NJError::ParseError(
                    "sequence data encountered before any FASTA header".into(),
                ));
            }
            current_seq.push_str(trimmed);
        }
    }

    // Flush last entry.
    if let Some(identifier) = current_name {
        if current_seq.is_empty() {
            return Err(NJError::ParseError(format!(
                "sequence '{identifier}' has no sequence data"
            )));
        }
        msa.push(SequenceObject { identifier, sequence: current_seq });
    }

    if msa.is_empty() {
        return Err(NJError::ParseError("input contains no FASTA sequences".into()));
    }

    // Alignment length check.
    let expected_len = msa[0].sequence.len();
    if expected_len == 0 {
        return Err(NJError::ParseError("sequences must not be empty".into()));
    }
    for seq in &msa[1..] {
        if seq.sequence.len() != expected_len {
            return Err(NJError::ParseError(format!(
                "sequence '{}' has length {}, expected {expected_len}",
                seq.identifier,
                seq.sequence.len()
            )));
        }
    }

    Ok(msa)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_fasta_simple() {
        let input = ">SeqA\nACGT\n>SeqB\nACGA\n";
        let seqs = parse_fasta(input).unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].identifier, "SeqA");
        assert_eq!(seqs[0].sequence, "ACGT");
        assert_eq!(seqs[1].identifier, "SeqB");
        assert_eq!(seqs[1].sequence, "ACGA");
    }

    #[test]
    fn test_parse_fasta_multiline_sequence() {
        let input = ">Seq1\nACGT\nACGT\n>Seq2\nAAAA\nAAAA\n";
        let seqs = parse_fasta(input).unwrap();
        assert_eq!(seqs[0].sequence, "ACGTACGT");
        assert_eq!(seqs[1].sequence, "AAAAAAAA");
    }

    #[test]
    fn test_parse_fasta_blank_lines_ignored() {
        let input = "\n>A\n\nACGT\n\n>B\nACGT\n\n";
        let seqs = parse_fasta(input).unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].sequence, "ACGT");
    }

    #[test]
    fn test_parse_fasta_header_trimmed() {
        let input = ">  My Sequence  \nACGT\n";
        let seqs = parse_fasta(input).unwrap();
        assert_eq!(seqs[0].identifier, "My Sequence");
    }

    #[test]
    fn test_parse_fasta_empty_input_errors() {
        assert!(parse_fasta("").is_err());
        assert!(parse_fasta("   \n\n  ").is_err());
    }

    #[test]
    fn test_parse_fasta_sequence_before_header_errors() {
        assert!(parse_fasta("ACGT\n>A\nACGT\n").is_err());
    }

    #[test]
    fn test_parse_fasta_header_with_no_sequence_errors() {
        assert!(parse_fasta(">A\n>B\nACGT\n").is_err());
        assert!(parse_fasta(">A\n").is_err());
    }

    #[test]
    fn test_parse_fasta_inconsistent_lengths_errors() {
        assert!(parse_fasta(">A\nACGT\n>B\nAC\n").is_err());
    }
}
