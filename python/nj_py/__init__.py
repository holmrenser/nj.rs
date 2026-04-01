from __future__ import annotations

from typing import Callable, Literal

from nj_py._nj_py import nj as _nj

SubstitutionModel = Literal["PDiff", "JukesCantor", "Kimura2P", "Poisson"]


def nj(
    msa: list[dict[str, str]],
    *,
    substitution_model: SubstitutionModel = "PDiff",
    n_bootstrap_samples: int = 0,
    on_progress: Callable[[int, int], None] | None = None,
) -> str:
    """Infer a Neighbor-Joining tree and return a Newick string.

    Args:
        msa: Aligned sequences as a list of ``{"identifier": ..., "sequence": ...}``
            dicts. All sequences must have the same length.
        substitution_model: Distance model. DNA: ``"PDiff"``, ``"JukesCantor"``,
            ``"Kimura2P"``. Protein: ``"PDiff"``, ``"Poisson"``. Alphabet is
            auto-detected; passing an incompatible model raises ``ValueError``.
        n_bootstrap_samples: Number of bootstrap replicates used to annotate
            internal nodes with support values. ``0`` disables bootstrapping.
        on_progress: Optional callable invoked as ``on_progress(completed, total)``
            after each bootstrap replicate. Never called when
            ``n_bootstrap_samples`` is 0.

    Returns:
        Newick string, e.g. ``"(A:0.1,B:0.2);"``

    Raises:
        ValueError: If the MSA is empty, a model is incompatible with the
            detected alphabet, or the NJ algorithm fails internally.

    Example::

        from tqdm import tqdm
        from nj_py import nj

        sequences = [
            {"identifier": "human", "sequence": "ACGTACGT"},
            {"identifier": "chimp", "sequence": "ACGCACGT"},
            {"identifier": "mouse", "sequence": "TCGTACGT"},
        ]

        with tqdm(total=100) as bar:
            tree = nj(
                sequences,
                substitution_model="JukesCantor",
                n_bootstrap_samples=100,
                on_progress=lambda current, total: bar.update(1),
            )
    """
    return _nj(
        {
            "msa": msa,
            "substitution_model": substitution_model,
            "n_bootstrap_samples": n_bootstrap_samples,
        },
        on_progress,
    )


__all__ = ["nj", "SubstitutionModel"]
