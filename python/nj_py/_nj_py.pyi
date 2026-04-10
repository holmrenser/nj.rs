from typing import Any, Callable

from nj_py import DistanceResult

def nj(
    py_config: dict[str, Any],
    on_progress: Callable[[int, int], None] | None = None,
) -> str: ...

def distance_matrix(py_config: dict[str, Any]) -> DistanceResult: ...

def average_distance(py_config: dict[str, Any]) -> float: ...
