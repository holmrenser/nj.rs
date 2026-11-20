//! Bootstrap support: replicate resampling, clade tallying and tree annotation.

use bitvec::prelude::{BitVec, Lsb0, bitvec};
use std::collections::HashMap;

use crate::alphabet::AlphabetEncoding;
use crate::config::MSA;
use crate::event::NJEvent;
use crate::models::{ModelCalculation, RateHet};
use crate::tree::{NameOrSupport, TreeNode};

/// Fills `out` with the leaf indices of all taxa in the subtree rooted at `node`.
///
/// Bits in `out` are set to `true` for each leaf encountered. The bit position
/// is looked up in `idx` by the leaf's name label. Returns `Err` if a leaf
/// has no name label (should not occur for well-formed NJ trees).
fn bitset_of(
    node: &TreeNode,
    idx: &HashMap<String, usize>,
    out: &mut BitVec<u8, Lsb0>,
) -> Result<(), String> {
    match &node.children {
        None => match &node.label {
            Some(NameOrSupport::Name(name)) => {
                let i = *idx
                    .get(name)
                    .ok_or_else(|| format!("leaf '{name}' not found in index map"))?;
                out.set(i, true);
                Ok(())
            }
            _ => Err("Leaf node without a name label".into()),
        },
        Some([l, r]) => {
            bitset_of(l, idx, out)?;
            bitset_of(r, idx, out)?;
            Ok(())
        }
    }
}

/// Recursively counts how many times each non-trivial clade appears in `tree`.
///
/// A clade is represented as a raw-byte encoding of a `BitVec` over the `n_taxa`
/// leaf indices. Only clades with `1 < size < n_taxa` (i.e. proper internal
/// clades) are counted. Each call increments the clade's entry in `counter` by 1.
/// Used by [`bootstrap_clade_counts`] to aggregate over bootstrap replicates.
pub(crate) fn count_clades(
    tree: &TreeNode,
    idx: &HashMap<String, usize>,
    n_taxa: usize,
    counter: &mut HashMap<Vec<u8>, usize>,
) -> Result<(), String> {
    if let Some([l, r]) = &tree.children {
        let mut bv = bitvec![u8, Lsb0; 0; n_taxa];
        bitset_of(tree, idx, &mut bv)?;

        let n = bv.count_ones();
        if n > 1 && n < n_taxa {
            counter
                .entry(bv.as_raw_slice().to_vec())
                .and_modify(|c| *c += 1)
                .or_insert(1);
        }

        count_clades(l, idx, n_taxa, counter)?;
        count_clades(r, idx, n_taxa, counter)?;
    }
    Ok(())
}

/// Runs a single bootstrap replicate end-to-end: resample alignment columns,
/// compute the distance matrix with model `M`, run Neighbor-Joining, and tally
/// the replicate's clades into a fresh counter.
///
/// Shared by the sequential and parallel bootstrap paths. Every fallible step is
/// propagated via `?` so a failed replicate becomes a returned `Err` rather than
/// a panic (which would otherwise poison a Rayon worker thread).
fn run_one_replicate<A, M>(
    msa: &MSA<A>,
    idx_map: &HashMap<String, usize>,
    n_taxa: usize,
    rates: RateHet,
) -> Result<HashMap<Vec<u8>, usize>, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    M: ModelCalculation<A> + Send + Sync,
{
    let tree = msa.bootstrap()?.into_dist::<M>(rates).neighbor_joining()?;
    let mut local = HashMap::new();
    count_clades(&tree, idx_map, n_taxa, &mut local)?;
    Ok(local)
}

/// Parallel bootstrap worker: runs all replicates on Rayon threads and sends
/// per-replicate clade maps over an MPSC channel. The main thread (caller)
/// merges results and fires `on_event` — keeping the callback on a single
/// thread with no `Sync` requirement.
///
/// Uses `std::thread::scope` so `msa` and `idx_map` can be borrowed into the
/// spawned thread without requiring `'static` lifetimes.
#[cfg(feature = "parallel")]
fn bootstrap_clade_counts_parallel<A, M>(
    msa: &MSA<A>,
    n_bootstrap_samples: usize,
    idx_map: &HashMap<String, usize>,
    n_taxa: usize,
    on_event: Option<&dyn Fn(NJEvent)>,
    num_threads: Option<usize>,
    rates: RateHet,
) -> Result<HashMap<Vec<u8>, usize>, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    M: ModelCalculation<A> + Send + Sync,
{
    use rayon::iter::{IntoParallelIterator, ParallelIterator};
    use std::sync::mpsc;

    let pool = crate::runner::build_thread_pool(num_threads)?;
    let (tx, rx) = mpsc::channel::<Result<HashMap<Vec<u8>, usize>, String>>();
    let mut counter: HashMap<Vec<u8>, usize> = HashMap::new();

    std::thread::scope(|scope| -> Result<(), String> {
        // Spawn the Rayon work onto a background OS thread so the main thread
        // can consume from `rx` concurrently, giving real per-replicate progress.
        // `pool.install` installs the thread pool as the active pool for the
        // duration of the closure, bounding parallelism to `num_threads` workers.
        scope.spawn(|| {
            pool.install(|| {
                (0..n_bootstrap_samples)
                    .into_par_iter()
                    .for_each_with(tx, |sender, _| {
                        let result = run_one_replicate::<A, M>(msa, idx_map, n_taxa, rates);
                        // Ignore send errors: only occurs if receiver was dropped,
                        // which cannot happen while we are in the recv loop below.
                        let _ = sender.send(result);
                    });
            });
        });

        // Main thread: receive results as they arrive, merge, and fire progress.
        for completed in 1..=n_bootstrap_samples {
            match rx.recv() {
                Ok(Ok(local)) => {
                    for (clade, count) in local {
                        *counter.entry(clade).or_insert(0) += count;
                    }
                }
                Ok(Err(e)) => return Err(e),
                Err(_) => return Err("bootstrap channel closed unexpectedly".into()),
            }
            if let Some(cb) = on_event {
                cb(NJEvent::BootstrapProgress {
                    completed,
                    total: n_bootstrap_samples,
                });
            }
        }
        Ok(())
    })?;

    Ok(counter)
}

/// Performs bootstrap sampling and counts clades across bootstrap trees.
///
/// When the `parallel` feature is enabled, replicates run on the Rayon thread
/// pool via [`bootstrap_clade_counts_parallel`]; the `on_event` callback is
/// still called after each replicate from the main thread. Without the feature
/// the loop is sequential and `on_event` fires in order.
pub(crate) fn bootstrap_clade_counts<A, M>(
    msa: &MSA<A>,
    n_bootstrap_samples: usize,
    on_event: Option<&dyn Fn(NJEvent)>,
    num_threads: Option<usize>,
    rates: RateHet,
) -> Result<Option<HashMap<Vec<u8>, usize>>, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    M: ModelCalculation<A> + Send + Sync,
{
    if n_bootstrap_samples == 0 {
        return Ok(None);
    }
    if let Some(cb) = on_event {
        cb(NJEvent::BootstrapStarted {
            total: n_bootstrap_samples,
        });
    }
    let idx_map: HashMap<String, usize> = msa.to_index_map();
    let n_taxa = msa.n_sequences;

    #[cfg(feature = "parallel")]
    let counter = bootstrap_clade_counts_parallel::<A, M>(
        msa,
        n_bootstrap_samples,
        &idx_map,
        n_taxa,
        on_event,
        num_threads,
        rates,
    )?;

    #[cfg(not(feature = "parallel"))]
    let counter = {
        let _ = num_threads;
        let mut c = HashMap::new();
        for i in 0..n_bootstrap_samples {
            let local = run_one_replicate::<A, M>(msa, &idx_map, n_taxa, rates)?;
            for (clade, count) in local {
                *c.entry(clade).or_insert(0) += count;
            }
            if let Some(cb) = on_event {
                cb(NJEvent::BootstrapProgress {
                    completed: i + 1,
                    total: n_bootstrap_samples,
                });
            }
        }
        c
    };

    Ok(Some(counter))
}

/// Annotates internal nodes with bootstrap support values from `counts`.
///
/// For each internal node, computes its clade `BitVec`, looks up the count in
/// `counts`, normalises it to a percentage (`count * 100 / n_bootstrap_samples`),
/// and assigns a [`NameOrSupport::Support`] label if a matching entry is found.
/// Nodes whose clade was never observed in bootstrap replicates receive no label.
pub(crate) fn add_bootstrap_to_tree(
    node: &mut TreeNode,
    idx: &HashMap<String, usize>,
    n_taxa: usize,
    counts: &HashMap<Vec<u8>, usize>,
    n_bootstrap_samples: usize,
) -> Result<(), String> {
    if node.children.is_some() {
        let mut bv = bitvec![u8, Lsb0; 0; n_taxa];
        bitset_of(node, idx, &mut bv)?;

        let n = bv.count_ones();
        if n > 1
            && n < n_taxa
            && let Some(c) = counts.get(&bv.as_raw_slice().to_vec())
        {
            let pct = c * 100 / n_bootstrap_samples;
            node.label = Some(NameOrSupport::Support(pct));
        }

        if let Some([l, r]) = &mut node.children {
            add_bootstrap_to_tree(l, idx, n_taxa, counts, n_bootstrap_samples)?;
            add_bootstrap_to_tree(r, idx, n_taxa, counts, n_bootstrap_samples)?;
        }
    }
    Ok(())
}
