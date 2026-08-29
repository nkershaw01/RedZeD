# redzed-tda

Fast Vietoris–Rips persistent homology, implemented in C++ with a NumPy interface.

This is the Python package for **RedZeD**, an algorithm to compute Vietoris-Rips
persistent homology. A key technique is *active enumeration*, allowing one to avoid enumerating 
almost all (n+1) birth simplices in a filtration.
RedZeD often gives a considerable improvement in both time and memory over 
existing implementations. Full details are in the paper:

> Chris Kapulkin and Nathan Kershaw.
> *RedZeD: Computing persistent homology by Reduction to Zero Differentials.*
> arXiv:2606.06310 (2026). <https://arxiv.org/abs/2606.06310>

The [main repository](https://github.com/nkershaw01/RedZeD) also holds the
reference Julia implementation.

If you encounter any errors or bugs, please email nkershaw@uwo.ca.

## Install

```bash
pip install redzed-tda
```

The distribution is `redzed-tda` and the import name is `redzed_tda`. (The name
`redzed` on PyPI belongs to an unrelated project.)

Wheels are provided for Linux, macOS and Windows, so no compiler is needed.

If no wheel matches your platform, or your Python version is newer than the last
release, pip falls back to building from source. That needs a C++17 compiler on
your machine. On Windows the compiler must be MSVC, not MinGW.

## Quick start

```python
import numpy as np
import redzed_tda as rz

points = np.random.default_rng(0).normal(size=(200, 3))
distances = rz.pairwise_distances(points)

h0, h1 = rz.rips(distances, maxdim=1)

print(h1)          # (k, 2) array of [birth, death] rows
```

## Input

`rips` takes a **distance matrix, not coordinates**. If your data is an array of 
points in R^n, either form works:

```python
# square (n, n) matrix
h0, h1 = rz.rips(rz.pairwise_distances(points), maxdim=1)

# condensed vector of n(n-1)/2 values, in scipy pdist order 
from scipy.spatial.distance import pdist
h0, h1 = rz.rips(pdist(points), maxdim=1)
```

pairwise_distances and pdist both take in an array of points and return the 
distance matrix. You may also pass redzed_tda.rips a square matrix constructed 
yourself. Only the lower triangle of a square input is read, so your matrix is
assumed symmetric with a zero diagonal. Passing a non-square array raises 
`ValueError`, but symmetry is not checked, so the user is responsible for 
ensuring the matrix is symmetric.

## Options

| Argument | Default | Meaning |
| --- | --- | --- |
| `distances` | — | Square `(n, n)` matrix or condensed `n(n-1)/2` vector. |
| `maxdim` | `1` | Highest homology dimension. `maxdim=1` returns H0 and H1. |
| `threshold` | `inf` | Ignore simplices with greater diameter. |
| `simplex_pairs` | `False` | Also return the birth and death simplex associated to every interval. |

## Output

A list of `maxdim + 1` arrays. Entry `d` has shape `(k, 2)` and holds the
`[birth, death]` pairs in dimension `d`:

```python
diagrams = rz.rips(distances, maxdim=2)
h0, h1, h2 = diagrams
```

With `simplex_pairs=True` you get `(diagrams, simplices)` instead, where
`simplices[d][i]` is a `(birth_vertices, death_vertices)` tuple corresponding to 
`diagrams[d][i]`. The death tuple is empty for an infinite interval.

## Building from source

```bash
git clone https://github.com/nkershaw01/RedZeD
cd RedZeD/python
pip install -e ".[test]"
pytest
```

Requires a C++17 compiler. On Windows that must be MSVC.

## Citing

If you use this software in published work, please cite the paper:

```bibtex
@misc{kapulkinkershaw2026redzed,
  title         = {RedZeD: Computing persistent homology by Reduction to Zero Differentials},
  author        = {Kapulkin, Chris and Kershaw, Nathan},
  year          = {2026},
  eprint        = {2606.06310},
  archivePrefix = {arXiv},
  primaryClass  = {cs.CG},
  doi           = {10.48550/arXiv.2606.06310},
  url           = {https://arxiv.org/abs/2606.06310}
}
```

## License

MIT. See [LICENSE](https://github.com/nkershaw01/RedZeD/blob/main/python/LICENSE).
