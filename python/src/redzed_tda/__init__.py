"""Fast Vietoris-Rips persistent homology.

    >>> import numpy as np, redzed_tda as rz
    >>> points = np.random.default_rng(0).normal(size=(50, 2))
    >>> diagrams = rz.rips(rz.pairwise_distances(points), maxdim=1)
    >>> diagrams[1].shape[1]
    2
"""

from __future__ import annotations

import math
from typing import List, Tuple, Union

import numpy as np

from . import _core

__all__ = ["rips", "pairwise_distances", "__version__"]

try:
    from importlib.metadata import PackageNotFoundError, version as _version

    __version__ = _version("redzed-tda")
except (ImportError, PackageNotFoundError):  # running from a source tree
    __version__ = "0.0.0.dev0"


def pairwise_distances(points: np.ndarray) -> np.ndarray:
    """Euclidean distance matrix for an (n_points, n_dimensions) array

    Provided so the package is useful without pulling in scipy 
    
    If you already have scipy, scipy.spatial.distance.pdist(points) 
    is equally acceptable input to rips
    """
    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2:
        raise ValueError(f"points must be a 2-D array, got a {points.ndim}-D array")
    # ||a - b||^2 = ||a||^2 - 2 a.b + ||b||^2, clipped because rounding can push
    # the diagonal a hair below zero
    sq = np.einsum("ij,ij->i", points, points)
    d2 = sq[:, None] - 2.0 * (points @ points.T) + sq[None, :]
    np.maximum(d2, 0.0, out=d2)
    d = np.sqrt(d2)
    np.fill_diagonal(d, 0.0)
    return d


def rips(
    distances: np.ndarray,
    maxdim: int = 1,
    threshold: float = math.inf,
    simplex_pairs: bool = False,
) -> Union[List[np.ndarray], Tuple[List[np.ndarray], List[list]]]:
    """Compute the Vietoris-Rips persistence diagrams of a distance matrix

    Parameters
    ----------
    distances:
        Either a square (n, n) distance matrix or a condensed 1-D vector of
        n(n-1)/2 values in scipy.spatial.distance.pdist order

        Only lower triangular entries read, user is responsible for ensuring the matrix is symmetric
    maxdim:
        Highest homology dimension 
    threshold:
        Radius cap
    simplex_pairs:
        Return the birth-death pairs of simplices as well

    Returns
    -------
    List of maxdim + 1 arrays, where entry d consists of [birth, death] intervals in dimension d.
    Deaths of inf mark features that never die 
    
    When simplex_pairs is true, returns
    (diagrams, simplices) instead, where simplices[d][i] is
    (birth_vertices, death_vertices) corresponding to the interval diagrams[d][i].
    For infinite intervals death_vertices will be empty. 
    """
    distances = np.ascontiguousarray(distances, dtype=np.float64)

    if not np.all(np.isfinite(distances)):
        raise ValueError("distance matrix contains NaN or infinite entries")
    if distances.size and distances.min() < 0.0:
        raise ValueError("distance matrix contains negative entries")
    if math.isnan(threshold):
        raise ValueError("threshold must be a number or inf, got nan")

    return _core.rips(
        distances,
        maxdim=int(maxdim),
        threshold=float(threshold),
        simplex_pairs=bool(simplex_pairs),
    )
