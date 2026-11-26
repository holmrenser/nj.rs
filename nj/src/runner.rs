//! Generic runners (distance matrix, average distance, full NJ) and the
//! runtime alphabet/model dispatch macro that instantiates them.

use std::collections::HashMap;

use crate::alphabet::AlphabetEncoding;
use crate::bootstrap::{add_bootstrap_to_tree, bootstrap_clade_counts};
use crate::config::{MSA, NJResult};
use crate::distance_matrix::DistanceResult;
use crate::event::{LogLevel, NJEvent};
use crate::models::{DistanceModel, RateHet};

/// Builds a Rayon thread pool with `num_threads` workers.
///
/// When `num_threads` is `None`, Rayon uses its default (one thread per logical CPU).
/// Returns `Err` if the pool cannot be constructed.
#[cfg(feature = "parallel")]
pub(crate) fn build_thread_pool(num_threads: Option<usize>) -> Result<rayon::ThreadPool, String> {
    let mut builder = rayon::ThreadPoolBuilder::new();
    if let Some(n) = num_threads {
        builder = builder.num_threads(n);
    }
    builder.build().map_err(|e| e.to_string())
}

/// Runs distance matrix computation with model `M` on alphabet `A`.
pub(crate) fn run_distance_matrix<A, D>(
    msa: MSA<A>,
    model: &D,
    num_threads: Option<usize>,
    rates: RateHet,
) -> Result<DistanceResult, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    D: DistanceModel<A>,
{
    #[cfg(feature = "parallel")]
    {
        let pool = build_thread_pool(num_threads)?;
        Ok(pool.install(|| msa.into_dist_with(model, rates)).into_result())
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = num_threads;
        Ok(msa.into_dist_with(model, rates).into_result())
    }
}

/// Runs average distance computation with model `M` on alphabet `A`.
pub(crate) fn run_average_distance<A, D>(
    msa: MSA<A>,
    model: &D,
    num_threads: Option<usize>,
    rates: RateHet,
) -> Result<f64, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    D: DistanceModel<A>,
{
    #[cfg(feature = "parallel")]
    {
        let pool = build_thread_pool(num_threads)?;
        Ok(pool.install(|| msa.into_dist_with(model, rates)).average())
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = num_threads;
        Ok(msa.into_dist_with(model, rates).average())
    }
}

/// Describes the pairs whose distance was capped, naming at most five of them.
fn saturation_message(names: &[String], capped: &[(usize, usize)], max_distance: f64) -> String {
    const SHOWN: usize = 5;
    let listed: Vec<String> = capped
        .iter()
        .take(SHOWN)
        .map(|&(i, j)| format!("'{}'/'{}'", names[j], names[i]))
        .collect();
    let more = capped.len().saturating_sub(SHOWN);
    let tail = if more > 0 {
        format!(" and {more} more")
    } else {
        String::new()
    };
    format!(
        "{} sequence pair(s) are saturated (no finite distance under the chosen model) and \
         were set to max_distance = {max_distance}: {}{tail}",
        capped.len(),
        listed.join(", ")
    )
}

/// Runs NJ with model `M` on alphabet `A` and returns an [`NJResult`].
///
/// If `n_bootstrap_samples > 0`, generates that many bootstrap replicates,
/// collects clade counts via [`bootstrap_clade_counts`], runs NJ on the
/// original distances, and annotates the tree before serialising to Newick.
/// When `include_distance_matrix` or `include_average_distance` is `true`,
/// the corresponding values are captured from the distance matrix before NJ
/// consumes it and included in the returned [`NJResult`].
pub(crate) fn run_nj<A, D>(
    msa: MSA<A>,
    model: &D,
    n_bootstrap_samples: usize,
    on_event: Option<&dyn Fn(NJEvent)>,
    num_threads: Option<usize>,
    include_distance_matrix: bool,
    include_average_distance: bool,
    rates: RateHet,
    max_distance: f64,
) -> Result<NJResult, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    D: DistanceModel<A>,
{
    let tally = bootstrap_clade_counts(
        model,
        &msa,
        n_bootstrap_samples,
        on_event,
        num_threads,
        rates,
        max_distance,
    )?;
    if let (Some(t), Some(cb)) = (&tally, on_event)
        && t.saturated_replicates > 0
    {
        cb(NJEvent::Log {
            level: LogLevel::Warning,
            message: format!(
                "{} of {} bootstrap replicates contained saturated sequence pairs; their \
                 distances were set to max_distance = {max_distance}",
                t.saturated_replicates, n_bootstrap_samples
            ),
        });
    }
    let clade_counts = tally.map(|t| t.counts);

    if let Some(cb) = on_event {
        cb(NJEvent::ComputingDistances);
    }

    #[cfg(feature = "parallel")]
    let mut dist = {
        let pool = build_thread_pool(num_threads)?;
        pool.install(|| msa.into_dist_with(model, rates))
    };
    #[cfg(not(feature = "parallel"))]
    let mut dist = msa.into_dist_with(model, rates);

    // Neighbor-joining cannot use infinite distances: cap them (and report) before
    // anything else reads the matrix, so the tree, the returned distance matrix and
    // the average all describe the same values.
    let capped = dist.cap_non_finite(max_distance);
    if let (false, Some(cb)) = (capped.is_empty(), on_event) {
        cb(NJEvent::Log {
            level: LogLevel::Warning,
            message: saturation_message(&dist.names, &capped, max_distance),
        });
    }

    let distance_matrix = if include_distance_matrix {
        Some(dist.to_result())
    } else {
        None
    };
    let average_distance = if include_average_distance {
        Some(dist.average())
    } else {
        None
    };

    if let Some(cb) = on_event {
        cb(NJEvent::RunningNJ);
    }

    let mut main_tree = dist.neighbor_joining()?;

    let newick = match clade_counts {
        Some(counts) => {
            if let Some(cb) = on_event {
                cb(NJEvent::AnnotatingBootstrap);
            }
            let main_idx_map: HashMap<String, usize> = msa.to_index_map();
            add_bootstrap_to_tree(
                &mut main_tree,
                &main_idx_map,
                msa.n_sequences,
                &counts,
                n_bootstrap_samples,
            )?;
            main_tree.to_newick()
        }
        None => main_tree.to_newick(),
    };
    Ok(NJResult {
        newick,
        distance_matrix,
        average_distance,
    })
}

/// Dispatches to a generic runner based on the runtime alphabet and substitution
/// model, enforcing the model–alphabet compatibility table in a single place.
///
/// `$matrix` is the `matrix` option (`Option<&str>`): required for, and only for, the `Custom`
/// model, whose rate matrix is parsed here. `$pairs` is consumed to build the typed [`MSA`];
/// `$run` is the generic runner (`run_nj`, `run_distance_matrix`, `run_average_distance`), called
/// as `$run(msa, &model, args…)` with the model as a value (the stateless models through
/// [`Static`](crate::models::Static)). Each runner returns `Result<_, String>`, mapped here to
/// [`NJError::AlgorithmFailure`]; unsupported (alphabet, model) pairs yield
/// [`NJError::IncompatibleModel`].
macro_rules! dispatch_run {
    ($alphabet:expr, $model:expr, $matrix:expr, $pairs:expr, $run:ident ( $($arg:expr),* $(,)? )) => {{
        use $crate::config::SubstitutionModel as Sm;
        use $crate::models::Static;
        let alphabet = $alphabet;
        let model = $model;
        let matrix: Option<&str> = $matrix;
        $crate::validate::validate_matrix(&model, matrix)?;
        let fail = $crate::error::NJError::AlgorithmFailure;
        match alphabet {
            $crate::alphabet::Alphabet::DNA => {
                let msa = $crate::config::MSA::<$crate::alphabet::DNA>::try_from_pairs($pairs).map_err(fail)?;
                match model {
                    Sm::PDiff => $run(msa, &Static::<$crate::models::PDiff>::new(), $($arg),*).map_err(fail),
                    Sm::JukesCantor => $run(msa, &Static::<$crate::models::JukesCantor>::new(), $($arg),*).map_err(fail),
                    Sm::Kimura2P => $run(msa, &Static::<$crate::models::Kimura2P>::new(), $($arg),*).map_err(fail),
                    Sm::TajimaNei => $run(msa, &Static::<$crate::models::TajimaNei>::new(), $($arg),*).map_err(fail),
                    Sm::Tamura => $run(msa, &Static::<$crate::models::Tamura>::new(), $($arg),*).map_err(fail),
                    Sm::Poisson | Sm::KimuraProtein | Sm::Custom => Err($crate::error::NJError::IncompatibleModel {
                        model,
                        alphabet: $crate::alphabet::Alphabet::DNA,
                    }),
                }
            }
            $crate::alphabet::Alphabet::Protein => {
                let msa = $crate::config::MSA::<$crate::alphabet::Protein>::try_from_pairs($pairs).map_err(fail)?;
                match model {
                    Sm::Poisson => $run(msa, &Static::<$crate::models::Poisson>::new(), $($arg),*).map_err(fail),
                    Sm::KimuraProtein => $run(msa, &Static::<$crate::models::KimuraProtein>::new(), $($arg),*).map_err(fail),
                    Sm::PDiff => $run(msa, &Static::<$crate::models::PDiff>::new(), $($arg),*).map_err(fail),
                    Sm::Custom => {
                        let text = matrix.expect("validate_matrix guarantees the matrix for Custom");
                        let custom = $crate::empirical::CustomProtein::from_paml(text)
                            .map_err(|reason| $crate::error::NJError::InvalidMatrix { reason })?;
                        $run(msa, &custom, $($arg),*).map_err(fail)
                    }
                    Sm::JukesCantor | Sm::Kimura2P | Sm::TajimaNei | Sm::Tamura => {
                        Err($crate::error::NJError::IncompatibleModel {
                            model,
                            alphabet: $crate::alphabet::Alphabet::Protein,
                        })
                    }
                }
            }
        }
    }};
}

pub(crate) use dispatch_run;
