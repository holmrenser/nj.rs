//! Command-line interface for the `nj` Neighbor-Joining tool.
//!
//! The CLI is gated behind the `cli` Cargo feature (which pulls in `clap`).
//! Build with `cargo build --features cli` or `cargo install nj --features cli`.
//!
//! # Usage
//!
//! ```text
//! nj [OPTIONS] [FASTA]
//!
//! Arguments:
//!   [FASTA]  MSA FASTA file to process (reads from stdin if omitted)
//!
//! Options:
//!   -b, --n-bootstrap-samples <N>    Number of bootstrap replicates [default: 100]
//!   -m, --substitution-model <MODEL> Substitution model [default: p-diff]
//!   -o, --output <FILE>              Write Newick output to file instead of stdout
//! ```

#[cfg(feature = "cli")]
mod cli {
    use ::nj::config::{DistConfig, NJConfig, SequenceObject};
    use ::nj::models::SubstitutionModel;
    use clap::Parser;
    use nj::nj;
    use std::fs;
    use std::io::Read;
    use std::path::PathBuf;

    /// Parsed command-line arguments.
    #[derive(Parser, Debug)]
    #[command(author, version, about)]
    pub struct Args {
        /// MSA FASTA file to process. Reads from stdin when omitted.
        #[arg(value_name = "FASTA")]
        pub input: Option<PathBuf>,

        /// Write Newick output to this file instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        pub output: Option<PathBuf>,

        /// Number of bootstrap replicates to generate (0 = no bootstrap).
        #[arg(short = 'b', long, default_value_t = 100)]
        pub n_bootstrap_samples: usize,

        /// Substitution model used to compute pairwise distances.
        #[arg(short = 'm', long, value_name = "MODEL", default_value = "p-diff")]
        pub substitution_model: SubstitutionModel,

        /// Output pairwise distance matrix as JSON instead of a Newick tree.
        #[arg(long, default_value_t = false)]
        pub distance_matrix: bool,

        /// Output the mean pairwise distance as a single number instead of a Newick tree.
        #[arg(long, default_value_t = false)]
        pub average_distance: bool,
    }

    /// Parses a FASTA-formatted string into a vector of [`SequenceObject`]s.
    ///
    /// - Whitespace at the start and end of each line is stripped.
    /// - Multi-line sequences are concatenated.
    /// - Empty lines are ignored.
    /// - Returns `Err` if a header has no sequence, a sequence appears before
    ///   any header, the file is empty, any sequence is empty, or sequences
    ///   differ in length.
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

    /// Reads the FASTA file from `args`, runs NJ, and writes the Newick
    /// output to stdout or a file.
    ///
    /// Argument parsing is handled by the caller; pass [`Args::parse()`] for
    /// production use or a manually constructed [`Args`] in tests.
    pub fn run(args: Args) -> Result<(), String> {
        use indicatif::{ProgressBar, ProgressStyle};

        if args.distance_matrix && args.average_distance {
            return Err("--distance-matrix and --average-distance are mutually exclusive".into());
        }

        let fasta = match &args.input {
            Some(path) => fs::read_to_string(path)
                .map_err(|e| format!("failed to read {}: {e}", path.display()))?,
            None => {
                let mut buf = String::new();
                std::io::stdin()
                    .read_to_string(&mut buf)
                    .map_err(|e| format!("failed to read stdin: {e}"))?;
                buf
            }
        };

        let msa = parse_fasta(&fasta)?;

        if args.distance_matrix {
            let conf = DistConfig {
                msa,
                substitution_model: args.substitution_model,
            };
            let result = ::nj::distance_matrix(conf)?;
            let json = serde_json::to_string_pretty(&result)
                .map_err(|e| format!("JSON serialization failed: {e}"))?;
            if let Some(path) = args.output {
                fs::write(&path, format!("{json}\n"))
                    .map_err(|e| format!("failed to write {}: {e}", path.display()))?;
            } else {
                println!("{json}");
            }
            return Ok(());
        }

        if args.average_distance {
            let conf = DistConfig {
                msa,
                substitution_model: args.substitution_model,
            };
            let avg = ::nj::average_distance(conf)?;
            if let Some(path) = args.output {
                fs::write(&path, format!("{avg}\n"))
                    .map_err(|e| format!("failed to write {}: {e}", path.display()))?;
            } else {
                println!("{avg}");
            }
            return Ok(());
        }

        let n_bootstrap = args.n_bootstrap_samples;
        let nj_conf = NJConfig {
            msa,
            n_bootstrap_samples: n_bootstrap,
            substitution_model: args.substitution_model,
        };

        let callback: Option<Box<dyn Fn(usize, usize)>> = if n_bootstrap > 0 {
            let pb = ProgressBar::new(n_bootstrap as u64);
            pb.set_style(
                ProgressStyle::with_template(
                    "[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} bootstrap",
                )
                .unwrap(),
            );
            Some(Box::new(move |current, _total| {
                pb.set_position(current as u64);
                if current == _total {
                    pb.finish_and_clear();
                }
            }))
        } else {
            None
        };

        let newick_tree = nj(nj_conf, callback)?;

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
        use clap::Parser;
        cli::run(cli::Args::parse())
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
    use super::cli::{run, Args};
    use ::nj::models::SubstitutionModel;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU32, Ordering};

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures")
            .join(name)
    }

    static COUNTER: AtomicU32 = AtomicU32::new(0);

    fn temp_fasta(content: &str) -> PathBuf {
        let n = COUNTER.fetch_add(1, Ordering::SeqCst);
        let path = std::env::temp_dir().join(format!("nj_test_{}_{}.fasta", std::process::id(), n));
        std::fs::write(&path, content).unwrap();
        path
    }

    fn base_args(input: PathBuf) -> Args {
        Args {
            input: Some(input),
            output: None,
            n_bootstrap_samples: 0,
            substitution_model: SubstitutionModel::PDiff,
            distance_matrix: false,
            average_distance: false,
        }
    }

    #[test]
    fn test_run_nj_succeeds() {
        assert!(run(base_args(fixture("simple_dna.fasta"))).is_ok());
    }

    #[test]
    fn test_run_distance_matrix_succeeds() {
        let args = Args {
            distance_matrix: true,
            ..base_args(fixture("simple_dna.fasta"))
        };
        assert!(run(args).is_ok());
    }

    #[test]
    fn test_run_average_distance_succeeds() {
        let args = Args {
            average_distance: true,
            ..base_args(fixture("simple_dna.fasta"))
        };
        assert!(run(args).is_ok());
    }

    #[test]
    fn test_run_writes_output_file() {
        let out = std::env::temp_dir().join(format!("nj_out_{}.nwk", std::process::id()));
        let args = Args {
            output: Some(out.clone()),
            ..base_args(fixture("simple_dna.fasta"))
        };
        run(args).expect("run failed");
        assert!(out.exists());
        let _ = std::fs::remove_file(out);
    }

    #[test]
    fn test_run_missing_file_is_error() {
        assert!(run(base_args(PathBuf::from("nonexistent.fasta"))).is_err());
    }

    #[test]
    fn test_run_empty_fasta_is_error() {
        let path = temp_fasta("");
        assert!(run(base_args(path)).is_err());
    }

    #[test]
    fn test_run_sequence_before_header_is_error() {
        let path = temp_fasta("ACGT\n>s\nAC\n");
        assert!(run(base_args(path)).is_err());
    }

    #[test]
    fn test_run_header_with_no_sequence_is_error() {
        let path = temp_fasta(">only_header\n");
        assert!(run(base_args(path)).is_err());
    }

    #[test]
    fn test_run_inconsistent_lengths_is_error() {
        let path = temp_fasta(">s1\nACGT\n>s2\nAC\n");
        assert!(run(base_args(path)).is_err());
    }

    // --- Conflicting output flags ---

    #[test]
    fn test_run_conflicting_output_flags_is_error() {
        let args = Args {
            distance_matrix: true,
            average_distance: true,
            ..base_args(fixture("simple_dna.fasta"))
        };
        assert!(run(args).is_err());
    }

    // --- Model–alphabet incompatibility ---

    #[test]
    fn test_run_dna_with_protein_model_is_error() {
        let args = Args {
            substitution_model: SubstitutionModel::Poisson,
            ..base_args(fixture("simple_dna.fasta"))
        };
        assert!(run(args).is_err());
    }

    #[test]
    fn test_run_protein_with_dna_only_model_is_error() {
        let args = Args {
            substitution_model: SubstitutionModel::JukesCantor,
            ..base_args(fixture("simple_protein.fasta"))
        };
        assert!(run(args).is_err());
    }

    // --- Non-default model success paths ---

    #[test]
    fn test_run_dna_jukes_cantor_succeeds() {
        let args = Args {
            substitution_model: SubstitutionModel::JukesCantor,
            ..base_args(fixture("simple_dna.fasta"))
        };
        assert!(run(args).is_ok());
    }

    #[test]
    fn test_run_dna_kimura2p_succeeds() {
        let args = Args {
            substitution_model: SubstitutionModel::Kimura2P,
            ..base_args(fixture("simple_dna.fasta"))
        };
        assert!(run(args).is_ok());
    }

    #[test]
    fn test_run_protein_poisson_succeeds() {
        let args = Args {
            substitution_model: SubstitutionModel::Poisson,
            ..base_args(fixture("simple_protein.fasta"))
        };
        assert!(run(args).is_ok());
    }

    #[test]
    fn test_run_protein_pdiff_succeeds() {
        assert!(run(base_args(fixture("simple_protein.fasta"))).is_ok());
    }
}
