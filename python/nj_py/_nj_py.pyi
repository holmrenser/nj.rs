from typing import Any, Callable

from nj_py import DistanceResult, NJEvent, NJResult

def nj(
    py_config: dict[str, Any],
    on_event: Callable[[NJEvent], None] | None = None,
) -> NJResult: ...

def distance_matrix(py_config: dict[str, Any]) -> DistanceResult: ...

def average_distance(py_config: dict[str, Any]) -> float: ...
