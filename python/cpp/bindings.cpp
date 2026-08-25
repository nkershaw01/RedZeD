// Python bindings for redzed

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "redzed.cpp"

namespace py = pybind11;

typedef py::array_t<double, py::array::c_style | py::array::forcecast> DoubleArray;

// ------------------------------------------------------------
// Input conversion
// ------------------------------------------------------------

// Accepts either a square (n, n) matrix or a condensed 1-D vector in scipy's
// pdist order
static DistanceMatrix matrix_from_array(const DoubleArray& arr) {
    DistanceMatrix D;

    if (arr.ndim() == 2) {
        if (arr.shape(0) != arr.shape(1))
            throw std::invalid_argument("a 2-D distance matrix must be square, got shape (" +
                                        std::to_string(arr.shape(0)) + ", " +
                                        std::to_string(arr.shape(1)) + ")");
        int nv = (int)arr.shape(0);
        D.resize_lower(nv);
        auto a = arr.unchecked<2>();
        // only the lower triangle is read, so an asymmetric input is not an error here
        for (int j = 1; j < nv; ++j) {
            value_t* dst = D.mat.data() + DistanceMatrix::offset(0, j);
            for (int i = 0; i < j; ++i) dst[i] = (value_t)a(j, i);
        }
    } else if (arr.ndim() == 1) {
        int64_t m = (int64_t)arr.shape(0);
        // invert m = nv(nv-1)/2
        int64_t nv = (int64_t)((1.0 + std::sqrt(1.0 + 8.0 * (double)m)) / 2.0 + 0.5);
        if (nv < 1 || nv * (nv - 1) / 2 != m)
            throw std::invalid_argument("a condensed distance vector must hold n(n-1)/2 values, "
                                        "and " + std::to_string(m) + " is not such a length");
        D.resize_lower((int)nv);
        auto a = arr.unchecked<1>();
        int64_t k = 0;
        for (int64_t i = 0; i < nv; ++i)
            for (int64_t j = i + 1; j < nv; ++j)
                D.mat[DistanceMatrix::offset((int)i, (int)j)] = (value_t)a((py::ssize_t)k++);
    } else {
        throw std::invalid_argument("expected a 2-D distance matrix or a 1-D condensed "
                                    "distance vector, got a " + std::to_string(arr.ndim()) +
                                    "-D array");
    }
    return D;
}

// ------------------------------------------------------------
// Output conversion
// ------------------------------------------------------------

// one array of [birth, death] intervals per dimension
static py::array_t<double> pairs_to_array(const std::vector<std::pair<value_t, value_t>>& v) {
    py::array_t<double> out({(py::ssize_t)v.size(), (py::ssize_t)2});
    auto r = out.mutable_unchecked<2>();
    for (size_t t = 0; t < v.size(); ++t) {
        r((py::ssize_t)t, 0) = (double)v[t].first;
        r((py::ssize_t)t, 1) = (double)v[t].second;
    }
    return out;
}

// (birth_vertices, death_vertices)
static py::list spx_to_list(const std::vector<SimplexPair>& v) {
    py::list out;
    for (const SimplexPair& sp : v)
        out.append(py::make_tuple(py::cast(sp.first), py::cast(sp.second)));
    return out;
}

// ------------------------------------------------------------
// Entry point
// ------------------------------------------------------------

static py::object rips_impl(const DoubleArray& distances, int maxdim, double threshold,
                            bool simplex_pairs) {
    if (maxdim < 0)
        throw std::invalid_argument("maxdim must be >= 0, got " + std::to_string(maxdim));

    DistanceMatrix D = matrix_from_array(distances);
    if (D.n < 1)
        throw std::invalid_argument("need at least one point");

    std::vector<std::vector<std::pair<value_t, value_t>>> pairs;
    std::vector<std::vector<SimplexPair>> spx;
    {
        // the solver never touches Python, so other threads can run meanwhile
        py::gil_scoped_release unlocked;
        pairs = redzed(D, maxdim, (value_t)threshold, simplex_pairs,
                       simplex_pairs ? &spx : nullptr);
    }

    py::list diagrams;
    for (int d = 0; d <= maxdim; ++d) diagrams.append(pairs_to_array(pairs[d]));
    if (!simplex_pairs) return diagrams;

    py::list simplices;
    for (int d = 0; d <= maxdim; ++d) simplices.append(spx_to_list(spx[d]));
    return py::make_tuple(diagrams, simplices);
}

PYBIND11_MODULE(_core, m) {
    m.doc() = "Compiled core of redzed: Vietoris-Rips persistent homology.";

    m.def("rips", &rips_impl,
          py::arg("distances"),
          py::arg("maxdim") = 1,
          py::arg("threshold") = std::numeric_limits<double>::infinity(),
          py::arg("simplex_pairs") = false,
          "Compute Vietoris-Rips persistence from a distance matrix.");
}
