use nj::{FastaSequence, MSA, NJConfig, nj};
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
    }

    pub fn run() -> Result<(), String> {
        let args = Args::parse();

        let fasta = fs::read_to_string(&args.input)
            .map_err(|e| format!("failed to read {}: {e}", args.input.display()))?;

        let msa = super::parse_fasta(&fasta)?;
        if msa.is_empty() {
            return Err("input FASTA contains no sequences".into());
        }
        let expected_len = msa[0].sequence.len();
        if expected_len == 0 {
            return Err("input FASTA sequences are empty".into());
        }
        for (i, fs) in msa.iter().enumerate() {
            if fs.sequence.len() != expected_len {
                return Err(format!(
                    "sequence {} ({}) has length {}, expected {}",
                    i,
                    fs.header,
                    fs.sequence.len(),
                    expected_len
                ));
            }
        }
        let newick_tree = nj(NJConfig {
            msa,
            hide_internal: true,
        })?;

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

fn parse_fasta(input: &str) -> Result<MSA, String> {
    let mut msa = MSA::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = String::new();

    for line in input.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(name) = current_name.replace(rest.trim().to_string()) {
                msa.push(FastaSequence {
                    header: name,
                    sequence: current_seq,
                });
                current_seq = String::new();
            }
        } else {
            if current_name.is_none() {
                return Err("FASTA sequence encountered before any header".into());
            }
            current_seq.push_str(trimmed);
        }
    }

    if let Some(name) = current_name {
        msa.push(FastaSequence {
            header: name,
            sequence: current_seq,
        });
    }

    Ok(msa)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_empty_input() {
        let msa = parse_fasta("").expect("should parse empty input");
        assert!(msa.is_empty());
    }

    #[test]
    fn test_parse_single_sequence() {
        let input = ">seq1\nACGT\n";
        let msa = parse_fasta(input).expect("parse failed");
        let mut expected = MSA::new();
        expected.push(FastaSequence {
            header: "seq1".into(),
            sequence: "ACGT".into(),
        });
        assert_eq!(msa.len(), expected.len());
        assert_eq!(msa[0].header, expected[0].header);
        assert_eq!(msa[0].sequence, expected[0].sequence);
    }

    #[test]
    fn test_parse_multiple_sequences_and_multiline_sequence() {
        let input = ">s1\nAA\nCC\n>s2\nGG\n";
        let msa = parse_fasta(input).expect("parse failed");
        let mut expected = MSA::new();
        expected.push(FastaSequence {
            header: "s1".into(),
            sequence: "AACC".into(),
        });
        expected.push(FastaSequence {
            header: "s2".into(),
            sequence: "GG".into(),
        });
        assert_eq!(msa.len(), expected.len());
        for (a, b) in msa.iter().zip(expected.iter()) {
            assert_eq!(a.header, b.header);
            assert_eq!(a.sequence, b.sequence);
        }
    }

    #[test]
    fn test_parse_with_blank_lines_and_trimming() {
        let input = "\n  > name \n   ACTG   \n\n";
        let msa = parse_fasta(input).expect("parse failed");
        let mut expected = MSA::new();
        expected.push(FastaSequence {
            header: "name".into(),
            sequence: "ACTG".into(),
        });
        assert_eq!(msa.len(), expected.len());
        assert_eq!(msa[0].header, expected[0].header);
        assert_eq!(msa[0].sequence, expected[0].sequence);
    }

    #[test]
    fn test_sequence_before_header_is_error() {
        let input = "ACGT\n>s\nAC\n";
        assert!(parse_fasta(input).is_err());
    }

    #[test]
    fn test_header_with_no_sequence_produces_empty_sequence() {
        let input = ">only_header\n";
        let msa = parse_fasta(input).expect("parse failed");
        assert_eq!(msa.len(), 1);
        assert_eq!(msa[0].header, "only_header");
        assert_eq!(msa[0].sequence, "");
    }

    #[test]
    fn test_neighbor_joining_three_taxa() {
        let input = ">A\nACG\n>B\nATG\n>C\nA-G\n";
        let msa = parse_fasta(input).expect("parse failed");
        let newick = nj(NJConfig {
            msa,
            hide_internal: false,
        })
        .expect("NJ failed");
        // The expected Newick string may vary depending on implementation details.
        // Here we just check that it contains the correct taxa names.
        assert!(newick.contains("A"));
        assert!(newick.contains("B"));
        assert!(newick.contains("C"));
    }
}
