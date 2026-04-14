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
//!   -v, --verbose                    Print stage log messages to stderr
//! ```

#[cfg(feature = "cli")]
mod cli {
    use ::nj::alphabet::Alphabet;
    use ::nj::config::NJConfig;
    use ::nj::event::{LogLevel, NJEvent};
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

        /// Write output to this file instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        pub output: Option<PathBuf>,

        /// Number of bootstrap replicates to generate (0 = no bootstrap).
        #[arg(short = 'b', long, default_value_t = 100)]
        pub n_bootstrap_samples: usize,

        /// Substitution model used to compute pairwise distances.
        #[arg(short = 'm', long, value_name = "MODEL", default_value = "p-diff")]
        pub substitution_model: SubstitutionModel,

        /// Override alphabet auto-detection. Detected from sequences when omitted.
        #[arg(short = 'a', long, value_name = "ALPHABET")]
        pub alphabet: Option<Alphabet>,

        /// Include the pairwise distance matrix in the output (switches output to JSON).
        #[arg(long, default_value_t = false)]
        pub distance_matrix: bool,

        /// Include the mean pairwise distance in the output (switches output to JSON).
        #[arg(long, default_value_t = false)]
        pub average_distance: bool,

        /// Number of threads to use for parallel computation (default: all available).
        /// Only effective when built with the `parallel` feature.
        #[arg(short = 't', long, value_name = "N")]
        pub num_threads: Option<usize>,

        /// Print algorithm stage messages (alphabet detection, distance computation, etc.) to stderr.
        #[arg(short = 'v', long, default_value_t = false)]
        pub verbose: bool,
    }

    /// Reads the FASTA file from `args`, runs NJ, and writes output to stdout or a file.
    ///
    /// When neither `--distance-matrix` nor `--average-distance` is set, output is a
    /// plain Newick string. When either flag is set, output is a JSON object whose
    /// `newick` field is always present, and `distance_matrix` / `average_distance`
    /// fields appear when the corresponding flags are set. Both flags may be combined.
    ///
    /// Argument parsing is handled by the caller; pass [`Args::parse()`] for
    /// production use or a manually constructed [`Args`] in tests.
    pub fn run(args: Args) -> Result<(), String> {
        use indicatif::{ProgressBar, ProgressStyle};

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

        let msa = ::nj::parse_fasta(&fasta).map_err(|e| e.to_string())?;
        let n_bootstrap = args.n_bootstrap_samples;
        let nj_conf = NJConfig {
            msa,
            n_bootstrap_samples: n_bootstrap,
            substitution_model: args.substitution_model,
            alphabet: args.alphabet,
            num_threads: args.num_threads,
            return_distance_matrix: args.distance_matrix,
            return_average_distance: args.average_distance,
        };

        let pb = if n_bootstrap > 0 {
            let pb = ProgressBar::new(n_bootstrap as u64);
            pb.set_style(
                ProgressStyle::with_template(
                    "[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} bootstrap",
                )
                .unwrap(),
            );
            Some(pb)
        } else {
            None
        };
        let verbose = args.verbose;
        let on_event: Box<dyn Fn(NJEvent)> = Box::new(move |event| match event {
            NJEvent::BootstrapProgress { completed, total } => {
                if let Some(ref pb) = pb {
                    pb.set_position(completed as u64);
                    if completed == total {
                        pb.finish_and_clear();
                    }
                }
            }
            NJEvent::Log { level, message } => {
                let tag = match level {
                    LogLevel::Info => "info",
                    LogLevel::Warning => "warning",
                };
                eprintln!("[{tag}] {message}");
            }
            stage_event => {
                if verbose {
                    let msg = match stage_event {
                        NJEvent::MsaValidated {
                            n_sequences,
                            n_sites,
                        } => {
                            format!("MSA validated: {n_sequences} sequences, {n_sites} sites")
                        }
                        NJEvent::AlphabetDetected { alphabet } => {
                            format!("Detected alphabet: {alphabet:?}")
                        }
                        NJEvent::ComputingDistances => "Computing distance matrix".to_string(),
                        NJEvent::RunningNJ => "Running Neighbor Joining".to_string(),
                        NJEvent::BootstrapStarted { total } => {
                            format!("Running {total} bootstrap replicates")
                        }
                        NJEvent::AnnotatingBootstrap => "Annotating bootstrap support".to_string(),
                        _ => return,
                    };
                    eprintln!("[info] {msg}");
                }
            }
        });

        let result = nj(nj_conf, Some(on_event)).map_err(|e| e.to_string())?;

        let output = if args.distance_matrix || args.average_distance {
            serde_json::to_string_pretty(&result)
                .map_err(|e| format!("JSON serialization failed: {e}"))?
        } else {
            result.newick
        };

        if let Some(path) = args.output {
            fs::write(&path, format!("{output}\n"))
                .map_err(|e| format!("failed to write {}: {e}", path.display()))?;
        } else {
            println!("{output}");
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
    use super::cli::{Args, run};
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
            alphabet: None,
            distance_matrix: false,
            average_distance: false,
            num_threads: None,
            verbose: false,
        }
    }

    #[test]
    fn test_run_all_outputs_succeeds() {
        let args = Args {
            distance_matrix: true,
            average_distance: true,
            ..base_args(fixture("simple_dna.fasta"))
        };
        assert!(run(args).is_ok());
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
        // Sequence data before first header is a parse error.
        let path = temp_fasta("ACGT\n>s\nACGT\n");
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
