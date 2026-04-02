import type { NJConfig, SequenceObject } from "./types/lib_types";

export type { NJConfig, SequenceObject };

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
 * @param onProgress - Optional callback invoked as `(current, total)` after
 *   each bootstrap replicate. Not called when `config.n_bootstrap_samples`
 *   is `0`.
 * @returns A semicolon-terminated Newick string. Bootstrap support counts (if
 *   requested) appear as integer labels on internal nodes.
 * @throws {Error} If the MSA is empty, sequences have unequal or zero length,
 *   or an incompatible model–alphabet combination is provided.
 */
export function nj(
  config: NJConfig,
  onProgress?: (current: number, total: number) => void,
): string;
