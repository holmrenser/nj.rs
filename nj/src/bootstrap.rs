//! Bootstrap support: replicate resampling, clade tallying and tree annotation.

use bitvec::prelude::{BitVec, Lsb0, bitvec};
use std::collections::{HashMap, HashSet};

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

/// Returns the rooting-independent key of the split induced by `node`, or `None`
/// for a leaf or a trivial split (one side has fewer than two taxa).
///
/// NJ trees are unrooted, but are stored with an arbitrary root (the final
/// join), so the same split can appear as clade `S` in one tree and as its
/// complement in another. The key is the raw-byte encoding of whichever side
/// does **not** contain taxon 0, which is identical for both orientations.
fn split_key(
    node: &TreeNode,
    idx: &HashMap<String, usize>,
    n_taxa: usize,
) -> Result<Option<Vec<u8>>, String> {
    if node.children.is_none() {
        return Ok(None);
    }
    let mut bv = bitvec![u8, Lsb0; 0; n_taxa];
    bitset_of(node, idx, &mut bv)?;

    let n = bv.count_ones();
    if n <= 1 || n + 1 >= n_taxa {
        return Ok(None);
    }
    if bv[0] {
        // Complement bit by bit so the unused padding bits stay zero and keys
        // compare equal as raw bytes.
        let mut flipped = bitvec![u8, Lsb0; 0; n_taxa];
        for i in 0..n_taxa {
            flipped.set(i, !bv[i]);
        }
        bv = flipped;
    }
    Ok(Some(bv.as_raw_slice().to_vec()))
}

/// Recursively collects the distinct splits of `tree` into `splits`.
///
/// Splits are keyed by [`split_key`]. A set is used so a split is recorded once
/// per tree: the two children of the root describe the same split, and must not
/// be counted twice. Used by [`bootstrap_clade_counts`] to aggregate over
/// bootstrap replicates.
pub(crate) fn collect_splits(
    tree: &TreeNode,
    idx: &HashMap<String, usize>,
    n_taxa: usize,
    splits: &mut HashSet<Vec<u8>>,
) -> Result<(), String> {
    if let Some([l, r]) = &tree.children {
        if let Some(key) = split_key(tree, idx, n_taxa)? {
            splits.insert(key);
        }
        collect_splits(l, idx, n_taxa, splits)?;
        collect_splits(r, idx, n_taxa, splits)?;
    }
    Ok(())
}

type Counts = HashMap<Vec<u8>, usize>;

/// Split tallies from the bootstrap phase.
pub(crate) struct BootstrapTally {
    /// How many replicates contained each split.
    pub counts: Counts,
    /// Number of replicates in which at least one distance was not finite and was
    /// replaced by `max_distance`.
    pub saturated_replicates: usize,
}

/// Runs a single bootstrap replicate end-to-end: resample alignment columns,
/// compute the distance matrix with model `M`, cap non-finite distances at
/// `max_distance`, run Neighbor-Joining, and tally the replicate's clades into a
/// fresh counter. The returned flag says whether any distance had to be capped
/// (resampling can saturate a pair that the original data do not).
///
/// Shared by the sequential and parallel bootstrap paths. Every fallible step is
/// propagated via `?` so a failed replicate becomes a returned `Err` rather than
/// a panic (which would otherwise poison a Rayon worker thread).
fn run_one_replicate<A, M>(
    msa: &MSA<A>,
    idx_map: &HashMap<String, usize>,
    n_taxa: usize,
    rates: RateHet,
    max_distance: f64,
) -> Result<(Counts, bool), String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    M: ModelCalculation<A> + Send + Sync,
{
    let mut dist = msa.bootstrap()?.into_dist::<M>(rates);
    let saturated = !dist.cap_non_finite(max_distance).is_empty();
    let tree = dist.neighbor_joining()?;
    let mut splits = HashSet::new();
    collect_splits(&tree, idx_map, n_taxa, &mut splits)?;
    Ok((splits.into_iter().map(|k| (k, 1)).collect(), saturated))
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
    max_distance: f64,
) -> Result<BootstrapTally, String>
where
    A: AlphabetEncoding + Send + Sync,
    A::Symbol: Send + Sync,
    M: ModelCalculation<A> + Send + Sync,
{
    use rayon::iter::{IntoParallelIterator, ParallelIterator};
    use std::sync::mpsc;

    let pool = crate::runner::build_thread_pool(num_threads)?;
    let (tx, rx) = mpsc::channel::<Result<(Counts, bool), String>>();
    let mut counter: Counts = HashMap::new();
    let mut saturated_replicates = 0usize;

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
                        let result = run_one_replicate::<A, M>(msa, idx_map, n_taxa, rates, max_distance);
                        // Ignore send errors: only occurs if receiver was dropped,
                        // which cannot happen while we are in the recv loop below.
                        let _ = sender.send(result);
                    });
            });
        });

        // Main thread: receive results as they arrive, merge, and fire progress.
        for completed in 1..=n_bootstrap_samples {
            match rx.recv() {
                Ok(Ok((local, saturated))) => {
                    for (clade, count) in local {
                        *counter.entry(clade).or_insert(0) += count;
                    }
                    saturated_replicates += usize::from(saturated);
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

    Ok(BootstrapTally {
        counts: counter,
        saturated_replicates,
    })
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
    max_distance: f64,
) -> Result<Option<BootstrapTally>, String>
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
    let tally = bootstrap_clade_counts_parallel::<A, M>(
        msa,
        n_bootstrap_samples,
        &idx_map,
        n_taxa,
        on_event,
        num_threads,
        rates,
        max_distance,
    )?;

    #[cfg(not(feature = "parallel"))]
    let tally = {
        let _ = num_threads;
        let mut c = HashMap::new();
        let mut saturated_replicates = 0usize;
        for i in 0..n_bootstrap_samples {
            let (local, saturated) =
                run_one_replicate::<A, M>(msa, &idx_map, n_taxa, rates, max_distance)?;
            for (clade, count) in local {
                *c.entry(clade).or_insert(0) += count;
            }
            saturated_replicates += usize::from(saturated);
            if let Some(cb) = on_event {
                cb(NJEvent::BootstrapProgress {
                    completed: i + 1,
                    total: n_bootstrap_samples,
                });
            }
        }
        BootstrapTally {
            counts: c,
            saturated_replicates,
        }
    };

    Ok(Some(tally))
}

/// Annotates internal nodes with bootstrap support values from `counts`.
///
/// For each internal node, computes its [`split_key`], looks up the count in
/// `counts`, normalises it to a percentage (`count * 100 / n_bootstrap_samples`),
/// and assigns a [`NameOrSupport::Support`] label if a matching entry is found.
/// Nodes whose split was never observed in bootstrap replicates receive no label.
pub(crate) fn add_bootstrap_to_tree(
    node: &mut TreeNode,
    idx: &HashMap<String, usize>,
    n_taxa: usize,
    counts: &HashMap<Vec<u8>, usize>,
    n_bootstrap_samples: usize,
) -> Result<(), String> {
    if node.children.is_some() {
        if let Some(key) = split_key(node, idx, n_taxa)?
            && let Some(c) = counts.get(&key)
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

#[cfg(test)]
mod tests {
    use super::*;

    fn leaf(i: usize, name: &str) -> Box<TreeNode> {
        Box::new(TreeNode::leaf(i, name.to_string(), Some(1.0)))
    }

    fn node(l: Box<TreeNode>, r: Box<TreeNode>) -> Box<TreeNode> {
        Box::new(TreeNode::internal(99, Some([l, r]), Some(1.0), None))
    }

    fn index() -> HashMap<String, usize> {
        ["A", "B", "C", "D", "E"]
            .iter()
            .enumerate()
            .map(|(i, n)| (n.to_string(), i))
            .collect()
    }

    fn splits_of(tree: &TreeNode) -> HashSet<Vec<u8>> {
        let mut s = HashSet::new();
        collect_splits(tree, &index(), 5, &mut s).unwrap();
        s
    }

    /// The unrooted tree ((A,B),(C,D),E) written with two different roots.
    #[test]
    fn same_unrooted_tree_gives_same_splits_for_different_roots() {
        // Rooted on the E branch: (((A,B),(C,D)),E)
        let t1 = node(
            node(
                node(leaf(0, "A"), leaf(1, "B")),
                node(leaf(2, "C"), leaf(3, "D")),
            ),
            leaf(4, "E"),
        );
        // Rooted on the (C,D) branch: ((A,B),E) joined to (C,D)
        let t2 = node(
            node(node(leaf(0, "A"), leaf(1, "B")), leaf(4, "E")),
            node(leaf(2, "C"), leaf(3, "D")),
        );
        // Rooted inside (A,B): (A,(B,((C,D),E)))
        let t3 = node(
            leaf(0, "A"),
            node(
                leaf(1, "B"),
                node(node(leaf(2, "C"), leaf(3, "D")), leaf(4, "E")),
            ),
        );
        let s1 = splits_of(&t1);
        assert_eq!(s1, splits_of(&t2));
        assert_eq!(s1, splits_of(&t3));
    }

    /// A split must be counted once per tree even though both root children
    /// describe it.
    #[test]
    fn root_split_is_counted_once() {
        let t = node(
            node(node(leaf(0, "A"), leaf(1, "B")), leaf(4, "E")),
            node(leaf(2, "C"), leaf(3, "D")),
        );
        // Splits: {A,B} | {C,D} | and {C,D}=={A,B,E}. Nontrivial distinct: AB, CD.
        assert_eq!(splits_of(&t).len(), 2);
    }

    #[test]
    fn different_topologies_give_different_splits() {
        let t1 = node(
            node(leaf(0, "A"), leaf(1, "B")),
            node(leaf(2, "C"), node(leaf(3, "D"), leaf(4, "E"))),
        );
        let t2 = node(
            node(leaf(0, "A"), leaf(2, "C")),
            node(leaf(1, "B"), node(leaf(3, "D"), leaf(4, "E"))),
        );
        assert_ne!(splits_of(&t1), splits_of(&t2));
    }

    /// Support is looked up correctly regardless of the root used to annotate.
    #[test]
    fn support_annotation_is_rooting_independent() {
        let mut counts = HashMap::new();
        for s in splits_of(&node(
            node(node(leaf(0, "A"), leaf(1, "B")), leaf(4, "E")),
            node(leaf(2, "C"), leaf(3, "D")),
        )) {
            counts.insert(s, 10usize);
        }
        // Different root, same unrooted topology.
        let mut t = node(
            leaf(4, "E"),
            node(
                node(leaf(0, "A"), leaf(1, "B")),
                node(leaf(2, "C"), leaf(3, "D")),
            ),
        );
        add_bootstrap_to_tree(&mut t, &index(), 5, &counts, 10).unwrap();
        let nwk = t.to_newick();
        assert_eq!(nwk.matches(")100").count(), 2, "{nwk}");
    }
}
