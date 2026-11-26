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

        /// Number of threads for parallel bootstrap and distance computation
        /// (default: all available cores). The CLI is always built with
        /// threading enabled, so this flag is always effective.
        #[arg(short = 't', long, value_name = "N")]
        pub num_threads: Option<usize>,

        /// Gamma rate-heterogeneity shape parameter (alpha, must be > 0). Omit for
        /// uniform substitution rates. Has no effect on the p-diff model.
        #[arg(short = 'g', long, value_name = "ALPHA")]
        pub gamma_shape: Option<f64>,

        /// Proportion of invariant sites in [0, 1). Omit to disable the
        /// invariant-sites correction. Has no effect on the p-diff model.
        #[arg(short = 'i', long = "p-invar", value_name = "FRACTION")]
        pub p_invar: Option<f64>,

        /// Rate matrix file for `--substitution-model custom` (and only then): PAML layout, the
        /// lower triangle of the symmetric exchangeability matrix by rows, then the 20
        /// amino-acid frequencies in the order ARNDCQEGHILKMFPSTWYV.
        #[arg(long, value_name = "FILE")]
        pub matrix: Option<PathBuf>,

        /// GTR exchangeabilities for `--substitution-model gtr` (and only then): six comma-separated
        /// relative values in the order AC,AG,AT,CG,CT,GT. HKY85 is `1,K,1,1,K,1`; all equal is F81.
        #[arg(long, value_name = "AC,AG,AT,CG,CT,GT", value_delimiter = ',', num_args = 1..)]
        pub gtr_rates: Option<Vec<f64>>,

        /// Base frequencies A,C,G,T (comma-separated) for the gtr model. Default: pooled over
        /// the alignment.
        #[arg(long, value_name = "A,C,G,T", value_delimiter = ',', num_args = 1..)]
        pub base_frequencies: Option<Vec<f64>>,

        /// Distance (substitutions per site) given to pairs whose distance is not finite
        /// because the sequences are saturated for the model. Default: 10.
        #[arg(long, value_name = "D")]
        pub max_distance: Option<f64>,

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
        let matrix = match &args.matrix {
            Some(path) => Some(
                fs::read_to_string(path)
                    .map_err(|e| format!("failed to read {}: {e}", path.display()))?,
            ),
            None => None,
        };
        let n_bootstrap = args.n_bootstrap_samples;
        let nj_conf = NJConfig {
            msa,
            n_bootstrap_samples: n_bootstrap,
            substitution_model: args.substitution_model,
            alphabet: args.alphabet,
            num_threads: args.num_threads,
            return_distance_matrix: args.distance_matrix,
            return_average_distance: args.average_distance,
            gamma_shape: args.gamma_shape,
            p_invar: args.p_invar,
            matrix,
            gtr_rates: args.gtr_rates,
            base_frequencies: args.base_frequencies,
            max_distance: args.max_distance,
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
            gamma_shape: None,
            p_invar: None,
            matrix: None,
            gtr_rates: None,
            base_frequencies: None,
            max_distance: None,
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

    #[test]
    fn test_run_protein_builtin_matrices_succeed() {
        for model in [
            SubstitutionModel::Lg,
            SubstitutionModel::Wag,
            SubstitutionModel::Jtt,
            SubstitutionModel::JttDcmut,
            SubstitutionModel::Dayhoff,
            SubstitutionModel::DayhoffDcmut,
        ] {
            let args = Args {
                substitution_model: model.clone(),
                ..base_args(fixture("simple_protein.fasta"))
            };
            assert!(run(args).is_ok(), "{model:?}");
        }
    }

    #[test]
    fn test_run_dna_gtr_succeeds_with_and_without_frequencies() {
        let args = Args {
            substitution_model: SubstitutionModel::Gtr,
            gtr_rates: Some(vec![1.0, 4.0, 1.0, 1.0, 4.0, 1.0]),
            distance_matrix: true,
            ..base_args(fixture("simple_dna.fasta"))
        };
        assert!(run(args).is_ok());
        let args = Args {
            substitution_model: SubstitutionModel::Gtr,
            gtr_rates: Some(vec![1.0, 4.0, 1.0, 1.0, 4.0, 1.0]),
            base_frequencies: Some(vec![0.3, 0.2, 0.2, 0.3]),
            ..base_args(fixture("simple_dna.fasta"))
        };
        assert!(run(args).is_ok());
    }

    #[test]
    fn test_run_gtr_option_errors() {
        let missing = Args {
            substitution_model: SubstitutionModel::Gtr,
            ..base_args(fixture("simple_dna.fasta"))
        };
        assert!(run(missing).unwrap_err().contains("gtr_rates"));
        let wrong_length = Args {
            substitution_model: SubstitutionModel::Gtr,
            gtr_rates: Some(vec![1.0, 2.0]),
            ..base_args(fixture("simple_dna.fasta"))
        };
        assert!(run(wrong_length).unwrap_err().contains("6 exchangeabilities"));
        let stray = Args {
            gtr_rates: Some(vec![1.0; 6]),
            ..base_args(fixture("simple_dna.fasta"))
        };
        assert!(run(stray).unwrap_err().contains("only used with the 'Gtr' model"));
        let protein = Args {
            substitution_model: SubstitutionModel::Gtr,
            gtr_rates: Some(vec![1.0; 6]),
            ..base_args(fixture("simple_protein.fasta"))
        };
        assert!(run(protein).is_err());
    }

    fn matrix_file(content: &str) -> PathBuf {
        let n = COUNTER.fetch_add(1, Ordering::SeqCst);
        let path = std::env::temp_dir().join(format!("nj_test_{}_{}.dat", std::process::id(), n));
        std::fs::write(&path, content).unwrap();
        path
    }

    #[test]
    fn test_run_protein_custom_matrix_succeeds() {
        let matrix = matrix_file(include_str!("../tests/data/equal_rates_protein.dat"));
        let args = Args {
            substitution_model: SubstitutionModel::Custom,
            matrix: Some(matrix),
            distance_matrix: true,
            ..base_args(fixture("simple_protein.fasta"))
        };
        assert!(run(args).is_ok());
    }

    #[test]
    fn test_run_custom_without_matrix_is_error() {
        let args = Args {
            substitution_model: SubstitutionModel::Custom,
            ..base_args(fixture("simple_protein.fasta"))
        };
        let err = run(args).unwrap_err();
        assert!(err.contains("matrix"), "{err}");
    }

    #[test]
    fn test_run_matrix_with_another_model_is_error() {
        let matrix = matrix_file(include_str!("../tests/data/equal_rates_protein.dat"));
        let args = Args {
            matrix: Some(matrix),
            ..base_args(fixture("simple_protein.fasta"))
        };
        assert!(run(args).is_err());
    }

    #[test]
    fn test_run_missing_or_malformed_matrix_file_is_error() {
        let missing = Args {
            substitution_model: SubstitutionModel::Custom,
            matrix: Some(PathBuf::from("/nonexistent/matrix.dat")),
            ..base_args(fixture("simple_protein.fasta"))
        };
        assert!(run(missing).unwrap_err().contains("failed to read"));
        let malformed = Args {
            substitution_model: SubstitutionModel::Custom,
            matrix: Some(matrix_file("1 2 3")),
            ..base_args(fixture("simple_protein.fasta"))
        };
        assert!(run(malformed).unwrap_err().contains("Invalid substitution matrix"));
    }
}
