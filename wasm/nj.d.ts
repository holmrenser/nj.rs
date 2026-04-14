import type { DistConfig, DistanceResult, NJConfig, NJEvent, NJResult, SequenceObject } from "./types/lib_types";

export type { DistConfig, DistanceResult, NJConfig, NJEvent, NJResult, SequenceObject };

export type NJEventCallback = (event: NJEvent) => void;

/**
 * Infers a Neighbor-Joining phylogenetic tree from a multiple sequence alignment.
 *
 * The alphabet (DNA vs. protein) is auto-detected from the sequences.
 * `config.substitution_model` must be compatible with the detected alphabet:
 * - DNA: `"PDiff"`, `"JukesCantor"`, `"Kimura2P"`
 * - Protein: `"PDiff"`, `"Poisson"`
 *
 * @param config - Alignment and run settings. All sequences in `config.msa`
 *   must have the same non-zero length. Set `config.n_bootstrap_samples` to
 *   `0` to skip bootstrapping.
 * @param onEvent - Optional callback invoked with a tagged event object at each
 *   algorithm stage. Use `event.type` to distinguish events. Bootstrap progress
 *   events have type `"BootstrapProgress"` with `completed` and `total` fields.
 * @returns An `NJResult` object. `result.newick` is always present as a
 *   semicolon-terminated Newick string. Bootstrap support values (if requested)
 *   appear as integer percentage labels (0–100) on internal nodes.
 *   `result.distance_matrix` and `result.average_distance` are present when
 *   the corresponding `return_distance_matrix` / `return_average_distance`
 *   flags are set in `config`.
 * @throws {Error} If the MSA is empty, sequences have unequal or zero length,
 *   or an incompatible model–alphabet combination is provided.
 */
export function nj(config: NJConfig, onEvent?: NJEventCallback): NJResult;

/**
 * Computes pairwise distances and returns a full symmetric distance matrix.
 *
 * @param config - Alignment and model settings. No bootstrap field needed.
 * @returns An object with `names` (taxon name array) and `matrix` (n×n
 *   symmetric distance matrix; diagonal entries are `0`).
 * @throws {Error} If the MSA is empty, sequences have unequal or zero length,
 *   or an incompatible model–alphabet combination is provided.
 */
export function distance_matrix(config: DistConfig): DistanceResult;

/**
 * Computes the mean of all n*(n-1)/2 unique pairwise distances.
 *
 * @param config - Alignment and model settings. No bootstrap field needed.
 * @returns Mean pairwise distance. Returns `0` for fewer than 2 taxa.
 * @throws {Error} If the MSA is empty, sequences have unequal or zero length,
 *   or an incompatible model–alphabet combination is provided.
 */
export function average_distance(config: DistConfig): number;
