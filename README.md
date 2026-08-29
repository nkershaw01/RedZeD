# RedZeD

A fast algorithm for Vietoris-Rips persistent homology, based on Reduction to Zero Differentials and active enumeration.

The code is associated with the paper [*RedZeD: Computing persistent homology by Reduction to Zero Differentials*](https://arxiv.org/abs/2606.06310) by Chris Kapulkin and Nathan Kershaw. 
The key innovation is active enumeration, which is made possible by Reduction to Zero Differentials. 
Using active enumeration, one is able to avoid enumerating almost all (n+1)-birth simplices, which make up the majority of all simplices in the filtration.
Due to this, RedZeD is considerably faster and more memory efficient than existing algorithms. 

There are three main things in this repo:

| What | Where | Notes |
| --- | --- | --- |
| Julia implementation | `redzed.jl` | Reference implementation |
| C++ implementation | `python/cpp/redzed.cpp` | Benchmarked against Ripser in the paper |
| Python package | `python/` | Published as `redzed-tda` |

## Python package

Found in the `python/` folder. 

Install using:
```bash
pip install redzed-tda
```
Example usage:
```python
import numpy as np, redzed_tda as rz
D = rz.pairwise_distances(points)
h0, h1 = rz.rips(D, maxdim=1)
```

Further details can be found in [python/README.md](python/README.md).

## Julia implementation

The Julia version is contained in the `redzed.jl` file and the primary function is redzed(dist, n; threshold = Inf, simplex_pairs = false). 


The inputs to `redzed(dist, n; threshold = Inf, simplex_pairs = false)` are:

- `dist`: the distance matrix 
    
- `n`: the max homology dimension 
    
- `threshold`: the radius threshold (default = Inf)

- `simplex_pairs`: return the birth/death simplex pairs associated to intervals? (default = false)
    
The output is a dictionary of the form: dimension => persistence pairs of that dimension.

If `simplex_pairs==true` it returns a tuple of dicts `(intervals,pairs)`, where `pairs[d][i]` is the simplex pair associated to the interval `intervals[d][i]`.

## C++ implementation

The C++ version is in `python/cpp/redzed.cpp`. 
It is the version tested against Ripser in the paper, and is the backend of the Python package.

### Building

```bash
g++ -O3 -std=c++17 -DNDEBUG -o redzed python/cpp/redzed.cpp           # Linux, macOS
g++ -O3 -std=c++17 -DNDEBUG -o redzed python/cpp/redzed.cpp -lpsapi   # Windows
```

`-lpsapi` is only needed for the memory reporting in `--timing`. MSVC works as well.

### Command line
Mainly just used for benchmarking. 

```
redzed <matrix.bin> <dim> [--lower] [--timing] [--threshold R] [--barcode] [--simplex-pairs]
```

| Flag | Meaning |
| --- | --- |
| `--lower` | read the input as a float32 lower triangle instead of a float64 square |
| `--threshold R` | ignore simplices above diameter `R` (default `+inf`) |
| `--barcode` | print every interval as `BAR <dim> <birth> <death>` |
| `--simplex-pairs` | print every interval together with its birth and death simplices |
| `--timing` | print `TIME` and, on Windows, `PEAKRSS` to stderr |

### Input

The input is binary, in one of two layouts:

| Mode | Layout |
| --- | --- |
| default | an `int64` point count, then that many rows of `float64`, row by row |
| `--lower` | raw `float32`, `n(n-1)/2` values, the lower triangle row by row, no header |

For default, only the lower triangle is read, the user is responsible for ensuring the matrix is symmetric.

To write lower matrix in python:
```python
import numpy as np
lower = D[np.tril_indices(len(D), k=-1)].astype(np.float32)
lower.tofile("matrix.bin")          # then: redzed matrix.bin 1 --lower
```

### Output

By default there is one summary line per dimension:

```
H0 512 1 44.5877756154
```

for each dimension this gives: the number of intervals, the number of infinite intervals, and a checksum of intervals (used for quick comparisons).
Use `--barcode` for the intervals:

```
BAR 1 0.1896519959 0.1975149959
BAR 0 0 inf
```

and `--simplex-pairs` to also get the simplices:

```
SPX 1 0.1896519959 0.1975149959 101 153 | 42 69 153
```

which is an H1 class born at 0.1897 by the edge (101, 153) and dying at 0.1975 when the triangle
(42, 69, 153) is added. 
Infinite intervals have no death simplex.

### Using it from C++

```cpp
#define REDZED_NO_MAIN
#include "redzed.cpp"

DistanceMatrix D;
D.resize_lower(n);
for (int j = 1; j < n; ++j)
    for (int i = 0; i < j; ++i)
        D.mat[DistanceMatrix::offset(i, j)] = my_distance(i, j);

auto pairs = redzed(D, 1, VINF);    // pairs[d] holds the (birth, death) intervals in dimension d
```

`REDZED_NO_MAIN` drops the command line driver. Note that `redzed` moves from `D`, so build a
fresh matrix for each call.

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

MIT. See [LICENSE](LICENSE).

