use nj::{NJConfig, SequenceObject, models::SubstitutionModel, nj};
use std::fs;
use std::path::PathBuf;

#[cfg(feature = "cli")]
mod cli {
    use super::*;
    use clap::Parser;

    #[derive(Parser, Debug)]
    #[command(author, version, about)]
    pub struct Args {
        /// MSA FASTA file to process
        #[arg(value_name = "FASTA")]
        pub input: PathBuf,

        /// Write Newick output to this file instead of stdout
        #[arg(short, long, value_name = "FILE")]
        pub output: Option<PathBuf>,

        /// Number of bootstrap samples to generate
        #[arg(short = 'b', long, default_value_t = 100)]
        pub n_bootstrap_samples: usize,

        /// Substitution model to use
        #[arg(short = 'm', long, value_name = "MODEL", default_value = "p-diff")]
        pub substitution_model: SubstitutionModel,
    }

    /// Parses a FASTA-formatted string into an MSA.
    /// Returns an error if the input is malformed.
    pub fn parse_fasta(input: &str) -> Result<Vec<SequenceObject>, String> {
        let mut msa = Vec::<SequenceObject>::new();
        let mut current_name: Option<String> = None;
        let mut current_seq = String::new();

        for line in input.lines() {
            // Trim whitespace
            let trimmed = line.trim();
            // Skip empty lines
            if trimmed.is_empty() {
                continue;
            }
            // Header line
            if let Some(fasta_header) = trimmed.strip_prefix('>') {
                // Save previous sequence if present
                if let Some(identifier) = current_name.replace(fasta_header.trim().to_string()) {
                    if current_seq.is_empty() {
                        return Err(format!(
                            "Sequence with identifier '{}' has no sequence",
                            identifier
                        ));
                    }
                    msa.push(SequenceObject {
                        identifier,
                        sequence: current_seq,
                    });
                    // Reset current sequence
                    current_seq = String::new();
                }
            // Sequence line
            } else {
                if current_name.is_none() {
                    return Err("FASTA sequence encountered before any header".into());
                }
                // Append to current sequence
                current_seq.push_str(trimmed);
            }
        }
        // Push the last sequence if present.
        if let Some(identifier) = current_name {
            if current_seq.is_empty() {
                return Err(format!(
                    "Sequence with identifier '{}' has no sequence",
                    identifier
                ));
            }
            msa.push(SequenceObject {
                identifier,
                sequence: current_seq,
            });
        }
        // Validate sequences
        if msa.is_empty() {
            return Err("input FASTA contains no sequences".into());
        }
        let expected_len = msa[0].len();
        if expected_len == 0 {
            return Err("input FASTA sequences are empty".into());
        }
        for (i, fs) in msa.iter().enumerate() {
            if fs.sequence.len() != expected_len {
                return Err(format!(
                    "sequence {} ({}) has length {}, expected {}",
                    i,
                    fs.identifier,
                    fs.sequence.len(),
                    expected_len
                ));
            }
        }

        Ok(msa)
    }

    /// Runs the CLI application.
    pub fn run() -> Result<(), String> {
        let args = Args::parse();

        let fasta = fs::read_to_string(&args.input)
            .map_err(|e| format!("failed to read {}: {e}", args.input.display()))?;

        let msa = parse_fasta(&fasta)?;

        let nj_conf = NJConfig {
            msa,
            n_bootstrap_samples: args.n_bootstrap_samples,
            substitution_model: args.substitution_model,
        };

        let newick_tree = nj(nj_conf)?;

        if let Some(path) = args.output {
            fs::write(&path, format!("{newick_tree}\n"))
                .map_err(|e| format!("failed to write {}: {e}", path.display()))?;
        } else {
            println!("{newick_tree}");
        }
        Ok(())
    }
}

fn main() -> Result<(), String> {
    #[cfg(feature = "cli")]
    {
        cli::run()
    }
    #[cfg(not(feature = "cli"))]
    {
        println!("CLI not enabled. Rebuild with --features cli");
        Ok(())
    }
}

#[cfg(test)]
#[cfg(feature = "cli")]
mod main_tests {
    use super::cli::parse_fasta;
    use super::*;

    #[test]
    fn test_parse_basic_fasta() {
        let input = ">seq1\nACGT\n>seq2\nTGCA\n";
        let expected = vec![
            SequenceObject {
                identifier: "seq1".into(),
                sequence: "ACGT".into(),
            },
            SequenceObject {
                identifier: "seq2".into(),
                sequence: "TGCA".into(),
            },
        ];
        let msa = parse_fasta(input).expect("parse failed");
        assert_eq!(msa.len(), expected.len());
        for (a, b) in msa.into_iter().zip(expected.into_iter()) {
            assert_eq!(a.identifier, b.identifier);
            assert_eq!(a.sequence, b.sequence);
        }
    }

    #[test]
    fn test_parse_empty_input_is_error() {
        assert!(parse_fasta("").is_err());
    }

    #[test]
    fn test_parse_single_sequence() {
        let input = ">s1\nAA\nCC\n";
        let expected = vec![SequenceObject {
            identifier: "s1".into(),
            sequence: "AACC".into(),
        }];
        let msa = parse_fasta(input).expect("parse failed");
        assert_eq!(msa.len(), expected.len());
        assert_eq!(msa[0].identifier, expected[0].identifier);
        assert_eq!(msa[0].sequence, expected[0].sequence);
    }

    #[test]
    fn test_parse_multiple_sequences_and_multiline_sequence() {
        let input = ">s1\nAA\nCC\n>s2\nGG\nTT\n";
        let expected = vec![
            SequenceObject {
                identifier: "s1".into(),
                sequence: "AACC".into(),
            },
            SequenceObject {
                identifier: "s2".into(),
                sequence: "GGTT".into(),
            },
        ];
        let msa = parse_fasta(input).expect("parse failed");

        assert_eq!(msa.len(), expected.len());
        for (a, b) in msa.into_iter().zip(expected.into_iter()) {
            assert_eq!(a.identifier, b.identifier);
            assert_eq!(a.sequence, b.sequence);
        }
    }

    #[test]
    fn test_parse_with_blank_lines_and_trimming() {
        let input = "\n  >s1    \nAA   \nCC\n\n>s2  \nGG\nTT\n";
        let expected = vec![
            SequenceObject {
                identifier: "s1".into(),
                sequence: "AACC".into(),
            },
            SequenceObject {
                identifier: "s2".into(),
                sequence: "GGTT".into(),
            },
        ];
        let msa = parse_fasta(input).expect("parse failed");
        assert_eq!(msa.len(), expected.len());
        assert_eq!(msa[0].identifier, expected[0].identifier);
        assert_eq!(msa[0].sequence, expected[0].sequence);
    }

    #[test]
    fn test_sequence_before_header_is_error() {
        let input = "ACGT\n>s\nAC\n";
        assert!(parse_fasta(input).is_err());
    }

    #[test]
    fn test_header_with_no_sequence_is_error() {
        let input = ">only_header\n";
        assert!(parse_fasta(input).is_err())
    }

    #[test]
    fn test_inconsistent_sequence_lengths_are_error() {
        let input = ">s1\nACGT\n>s2\nAC\n";
        assert!(parse_fasta(input).is_err());
    }

    #[test]
    fn test_empty_sequence_is_error() {
        let input = ">s1\n\n>s2\nACGT\n";
        assert!(parse_fasta(input).is_err());
    }

    #[test]
    fn test_parse_fasta_with_whitespace() {
        let input = "   >seq1   \n  ACGT  \n>seq2\n TGCA \n ";
        let expected = vec![
            SequenceObject {
                identifier: "seq1".into(),
                sequence: "ACGT".into(),
            },
            SequenceObject {
                identifier: "seq2".into(),
                sequence: "TGCA".into(),
            },
        ];
        let msa = parse_fasta(input).expect("parse failed");
        assert_eq!(msa.len(), expected.len());
        for (a, b) in msa.into_iter().zip(expected.into_iter()) {
            assert_eq!(a.identifier, b.identifier);
            assert_eq!(a.sequence, b.sequence);
        }
    }

    #[test]
    fn test_parse_fasta_with_no_sequences_is_error() {
        let input = ">seq1\n>seq2\n";
        assert!(parse_fasta(input).is_err());
    }

    #[test]
    fn test_nj_single_taxon() {
        let input = ">A\nACGT\n";
        let msa = parse_fasta(input).expect("parse failed");
        let newick = nj(NJConfig {
            msa,
            n_bootstrap_samples: 1,
            substitution_model: SubstitutionModel::PDiff,
        })
        .expect("NJ failed");
        assert_eq!(newick, "A;");
    }

    #[test]
    fn test_nj_two_taxa() {
        let input = ">A\nACG\n>B\nATG\n";
        let msa = parse_fasta(input).expect("parse failed");
        let newick = nj(NJConfig {
            msa,
            n_bootstrap_samples: 1,
            substitution_model: SubstitutionModel::PDiff,
        })
        .expect("NJ failed");
        assert_eq!(newick, "(A:0.167,B:0.167);");
    }

    #[test]
    fn test_nj_three_taxa() {
        let input = ">A\nACG\n>B\nATG\n>C\nA-G\n";
        let msa = parse_fasta(input).expect("parse failed");
        let newick = nj(NJConfig {
            msa,
            n_bootstrap_samples: 1,
            substitution_model: SubstitutionModel::PDiff,
        })
        .expect("NJ failed");
        // The expected Newick string may vary depending on implementation details.
        // Here we just check that it contains the correct taxa names.
        assert!(newick.contains("A"));
        assert!(newick.contains("B"));
        assert!(newick.contains("C"));
    }

    #[test]
    fn test_nj_protein_sequences() {
        let input = ">Prot1\nACDEFGHIK\n>Prot2\nACDFFGHIK\n";
        let msa = parse_fasta(input).expect("parse failed");
        let newick = nj(NJConfig {
            msa,
            n_bootstrap_samples: 1,
            substitution_model: SubstitutionModel::PDiff,
        })
        .expect("NJ failed");
        assert!(newick.contains("Prot1"));
        assert!(newick.contains("Prot2"));
    }

    #[test]
    fn test_nj_empty_msa_is_error() {
        let msa = Vec::<SequenceObject>::new();
        let result = nj(NJConfig {
            msa,
            n_bootstrap_samples: 1,
            substitution_model: SubstitutionModel::PDiff,
        });
        assert!(result.is_err());
    }
}
