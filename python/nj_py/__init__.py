from __future__ import annotations

from typing import Callable, Literal, TypedDict

from nj_py._nj_py import average_distance as _average_distance
from nj_py._nj_py import distance_matrix as _distance_matrix
from nj_py._nj_py import nj as _nj

SubstitutionModel = Literal["PDiff", "JukesCantor", "Kimura2P", "Poisson"]

class SequenceObject(TypedDict):
    identifier: str
    sequence: str
class DistanceResult(TypedDict):
    names: list[str]
    matrix: list[list[float]]


def nj(
    msa: list[SequenceObject],
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


def distance_matrix(
    msa: list[SequenceObject],
    *,
    substitution_model: SubstitutionModel = "PDiff",
) -> DistanceResult:
    """Compute pairwise distances and return a full symmetric distance matrix.

    Args:
        msa: Aligned sequences as a list of ``{"identifier": ..., "sequence": ...}``
            dicts. All sequences must have the same length.
        substitution_model: Distance model. DNA: ``"PDiff"``, ``"JukesCantor"``,
            ``"Kimura2P"``. Protein: ``"PDiff"``, ``"Poisson"``. Alphabet is
            auto-detected; passing an incompatible model raises ``ValueError``.

    Returns:
        A dict with keys ``"names"`` (list of taxon names) and ``"matrix"``
        (n×n list of lists of floats, symmetric, diagonal zero).

    Raises:
        ValueError: If the MSA is empty, a model is incompatible with the
            detected alphabet, or sequences have unequal length.

    Example::

        from nj_py import distance_matrix

        sequences = [
            {"identifier": "human", "sequence": "ACGTACGT"},
            {"identifier": "chimp", "sequence": "ACGCACGT"},
            {"identifier": "mouse", "sequence": "TCGTACGT"},
        ]

        result = distance_matrix(sequences)
        # result["names"] == ["human", "chimp", "mouse"]
        # result["matrix"][0][1]  # distance between human and chimp
    """
    return _distance_matrix({"msa": msa, "substitution_model": substitution_model})


def average_distance(
    msa: list[SequenceObject],
    *,
    substitution_model: SubstitutionModel = "PDiff",
) -> float:
    """Compute the mean of all n*(n-1)/2 unique pairwise distances.

    Args:
        msa: Aligned sequences as a list of ``{"identifier": ..., "sequence": ...}``
            dicts. All sequences must have the same length.
        substitution_model: Distance model. DNA: ``"PDiff"``, ``"JukesCantor"``,
            ``"Kimura2P"``. Protein: ``"PDiff"``, ``"Poisson"``. Alphabet is
            auto-detected; passing an incompatible model raises ``ValueError``.

    Returns:
        Mean pairwise distance as a float. Returns ``0.0`` for fewer than 2 taxa.

    Raises:
        ValueError: If the MSA is empty, a model is incompatible with the
            detected alphabet, or sequences have unequal length.
    """
    return _average_distance({"msa": msa, "substitution_model": substitution_model})


__all__ = [
    "nj",
    "distance_matrix",
    "average_distance",
    "SubstitutionModel",
    "SequenceObject",
    "DistanceResult",
]
