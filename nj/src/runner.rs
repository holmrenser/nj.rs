//! Generic runners (distance matrix, average distance, full NJ) and the
//! runtime alphabet/model dispatch macro that instantiates them.

use std::collections::HashMap;

use crate::alphabet::AlphabetEncoding;
use crate::bootstrap::{add_bootstrap_to_tree, bootstrap_clade_counts};
use crate::config::{MSA, NJResult};
use crate::distance_matrix::DistanceResult;
use crate::event::NJEvent;
use crate::models::{ModelCalculation, RateHet};

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
pub(crate) fn run_distance_matrix<A, M>(
    msa: MSA<A>,
    num_threads: Option<usize>,
    rates: RateHet,
) -> Result<DistanceResult, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    M: ModelCalculation<A> + Send + Sync,
{
    #[cfg(feature = "parallel")]
    {
        let pool = build_thread_pool(num_threads)?;
        Ok(pool.install(|| msa.into_dist::<M>(rates)).into_result())
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = num_threads;
        Ok(msa.into_dist::<M>(rates).into_result())
    }
}

/// Runs average distance computation with model `M` on alphabet `A`.
pub(crate) fn run_average_distance<A, M>(
    msa: MSA<A>,
    num_threads: Option<usize>,
    rates: RateHet,
) -> Result<f64, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    M: ModelCalculation<A> + Send + Sync,
{
    #[cfg(feature = "parallel")]
    {
        let pool = build_thread_pool(num_threads)?;
        Ok(pool.install(|| msa.into_dist::<M>(rates)).average())
    }
    #[cfg(not(feature = "parallel"))]
    {
        let _ = num_threads;
        Ok(msa.into_dist::<M>(rates).average())
    }
}

/// Runs NJ with model `M` on alphabet `A` and returns an [`NJResult`].
///
/// If `n_bootstrap_samples > 0`, generates that many bootstrap replicates,
/// collects clade counts via [`bootstrap_clade_counts`], runs NJ on the
/// original distances, and annotates the tree before serialising to Newick.
/// When `include_distance_matrix` or `include_average_distance` is `true`,
/// the corresponding values are captured from the distance matrix before NJ
/// consumes it and included in the returned [`NJResult`].
pub(crate) fn run_nj<A, M>(
    msa: MSA<A>,
    n_bootstrap_samples: usize,
    on_event: Option<&dyn Fn(NJEvent)>,
    num_threads: Option<usize>,
    include_distance_matrix: bool,
    include_average_distance: bool,
    rates: RateHet,
) -> Result<NJResult, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    M: ModelCalculation<A> + Send + Sync,
{
    let clade_counts =
        bootstrap_clade_counts::<A, M>(&msa, n_bootstrap_samples, on_event, num_threads, rates)?;

    if let Some(cb) = on_event {
        cb(NJEvent::ComputingDistances);
    }

    #[cfg(feature = "parallel")]
    let dist = {
        let pool = build_thread_pool(num_threads)?;
        pool.install(|| msa.into_dist::<M>(rates))
    };
    #[cfg(not(feature = "parallel"))]
    let dist = msa.into_dist::<M>(rates);

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
/// `$pairs` is consumed to build the typed [`MSA`]; `$run` is the generic runner
/// (`run_nj`, `run_distance_matrix`, `run_average_distance`) and its trailing
/// arguments follow the `msa` it receives. Each runner returns `Result<_, String>`,
/// mapped here to [`NJError::AlgorithmFailure`]; unsupported (alphabet, model)
/// pairs yield [`NJError::IncompatibleModel`].
macro_rules! dispatch_run {
    ($alphabet:expr, $model:expr, $pairs:expr, $run:ident ( $($arg:expr),* $(,)? )) => {{
        let alphabet = $alphabet;
        let model = $model;
        match alphabet {
            $crate::alphabet::Alphabet::DNA => {
                let msa = $crate::config::MSA::<$crate::alphabet::DNA>::try_from_pairs($pairs).map_err($crate::error::NJError::AlgorithmFailure)?;
                match model {
                    $crate::config::SubstitutionModel::PDiff => {
                        $run::<$crate::alphabet::DNA, $crate::models::PDiff>(msa, $($arg),*).map_err($crate::error::NJError::AlgorithmFailure)
                    }
                    $crate::config::SubstitutionModel::JukesCantor => {
                        $run::<$crate::alphabet::DNA, $crate::models::JukesCantor>(msa, $($arg),*).map_err($crate::error::NJError::AlgorithmFailure)
                    }
                    $crate::config::SubstitutionModel::Kimura2P => {
                        $run::<$crate::alphabet::DNA, $crate::models::Kimura2P>(msa, $($arg),*).map_err($crate::error::NJError::AlgorithmFailure)
                    }
                    $crate::config::SubstitutionModel::TajimaNei => {
                        $run::<$crate::alphabet::DNA, $crate::models::TajimaNei>(msa, $($arg),*).map_err($crate::error::NJError::AlgorithmFailure)
                    }
                    $crate::config::SubstitutionModel::Tamura => {
                        $run::<$crate::alphabet::DNA, $crate::models::Tamura>(msa, $($arg),*).map_err($crate::error::NJError::AlgorithmFailure)
                    }
                    $crate::config::SubstitutionModel::Poisson | $crate::config::SubstitutionModel::KimuraProtein => {
                        Err($crate::error::NJError::IncompatibleModel {
                            model,
                            alphabet: $crate::alphabet::Alphabet::DNA,
                        })
                    }
                }
            }
            $crate::alphabet::Alphabet::Protein => {
                let msa =
                    $crate::config::MSA::<$crate::alphabet::Protein>::try_from_pairs($pairs).map_err($crate::error::NJError::AlgorithmFailure)?;
                match model {
                    $crate::config::SubstitutionModel::Poisson => {
                        $run::<$crate::alphabet::Protein, $crate::models::Poisson>(msa, $($arg),*).map_err($crate::error::NJError::AlgorithmFailure)
                    }
                    $crate::config::SubstitutionModel::KimuraProtein => {
                        $run::<$crate::alphabet::Protein, $crate::models::KimuraProtein>(msa, $($arg),*)
                            .map_err($crate::error::NJError::AlgorithmFailure)
                    }
                    $crate::config::SubstitutionModel::PDiff => {
                        $run::<$crate::alphabet::Protein, $crate::models::PDiff>(msa, $($arg),*).map_err($crate::error::NJError::AlgorithmFailure)
                    }
                    $crate::config::SubstitutionModel::JukesCantor
                    | $crate::config::SubstitutionModel::Kimura2P
                    | $crate::config::SubstitutionModel::TajimaNei
                    | $crate::config::SubstitutionModel::Tamura => Err($crate::error::NJError::IncompatibleModel {
                        model,
                        alphabet: $crate::alphabet::Alphabet::Protein,
                    }),
                }
            }
        }
    }};
}

pub(crate) use dispatch_run;
