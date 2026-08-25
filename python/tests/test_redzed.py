import math

import numpy as np
import pytest

import redzed_tda as rz


def circle(n, radius=1.0):
    t = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return np.column_stack([radius * np.cos(t), radius * np.sin(t)])


def test_two_points_have_one_finite_h0_interval():
    d = np.array([[0.0, 2.0], [2.0, 0.0]])
    h0, _ = rz.rips(d, maxdim=1)
    finite = h0[np.isfinite(h0[:, 1])]
    assert finite.shape == (1, 2)
    assert finite[0, 0] == 0.0
    assert finite[0, 1] == pytest.approx(2.0)


def test_circle_has_one_long_h1_feature():
    """A sampled circle is the canonical loop: exactly one prominent H1 class."""
    d = rz.pairwise_distances(circle(40))
    h0, h1 = rz.rips(d, maxdim=1)

    # 40 points, one component survives forever, the other 39 die
    assert np.isinf(h0[:, 1]).sum() == 1
    assert h0.shape[0] == 40

    lifetimes = h1[:, 1] - h1[:, 0]
    assert h1.shape[0] >= 1
    longest = np.sort(lifetimes)[-1]
    runner_up = np.sort(lifetimes)[-2] if len(lifetimes) > 1 else 0.0
    assert longest > 10 * max(runner_up, 1e-12)


def test_condensed_and_square_inputs_agree():
    points = np.random.default_rng(0).normal(size=(30, 3))
    square = rz.pairwise_distances(points)
    condensed = square[np.triu_indices(30, k=1)]

    from_square = rz.rips(square, maxdim=1)
    from_condensed = rz.rips(condensed, maxdim=1)

    for a, b in zip(from_square, from_condensed):
        np.testing.assert_allclose(np.sort(a, axis=0), np.sort(b, axis=0), rtol=1e-5)


def test_threshold_reduces_or_preserves_interval_count():
    d = rz.pairwise_distances(circle(30))
    full = rz.rips(d, maxdim=1)
    capped = rz.rips(d, maxdim=1, threshold=0.5)
    assert capped[1].shape[0] <= full[1].shape[0]


def test_maxdim_controls_number_of_diagrams():
    d = rz.pairwise_distances(circle(15))
    assert len(rz.rips(d, maxdim=0)) == 1
    assert len(rz.rips(d, maxdim=1)) == 2
    assert len(rz.rips(d, maxdim=2)) == 3


@pytest.mark.parametrize("seed", range(5))
def test_maxdim_zero_agrees_with_higher_maxdim(seed):
    """maxdim=0 takes a separate branch in the solver, so check it against maxdim=1.

    This is the case that used to read R[1] out of bounds.
    """
    points = np.random.default_rng(seed).normal(size=(60, 3))
    d = rz.pairwise_distances(points)
    only_h0 = rz.rips(d, maxdim=0)[0]
    with_h1 = rz.rips(d, maxdim=1)[0]
    np.testing.assert_array_equal(np.sort(only_h0, axis=0), np.sort(with_h1, axis=0))


def test_maxdim_zero_gives_one_component_per_cluster():
    """Three well separated clusters, cut below the gap: three infinite H0 classes."""
    rng = np.random.default_rng(3)
    points = np.vstack([rng.normal(loc, 0.05, size=(15, 2)) for loc in ([0, 0], [10, 0], [0, 10])])
    (h0,) = rz.rips(rz.pairwise_distances(points), maxdim=0, threshold=1.0)
    assert np.isinf(h0[:, 1]).sum() == 3


def test_simplex_pairs_line_up_with_diagrams():
    d = rz.pairwise_distances(circle(20))
    diagrams, simplices = rz.rips(d, maxdim=1, simplex_pairs=True)
    assert len(diagrams) == len(simplices)
    for dim, (diag, spx) in enumerate(zip(diagrams, simplices)):
        assert len(spx) == diag.shape[0]
        for (birth_verts, death_verts), (birth, death) in zip(spx, diag):
            assert len(birth_verts) == dim + 1
            # an infinite interval has no death simplex
            assert (len(death_verts) == 0) == bool(math.isinf(death))


def test_births_never_exceed_deaths():
    points = np.random.default_rng(7).normal(size=(40, 2))
    for diagram in rz.rips(rz.pairwise_distances(points), maxdim=1):
        assert np.all(diagram[:, 0] <= diagram[:, 1])


@pytest.mark.parametrize(
    "bad, message",
    [
        (np.array([[0.0, 1.0]]), "square"),
        (np.array([0.0, 1.0]), "n\\(n-1\\)/2"),
        (np.array([[0.0, np.nan], [np.nan, 0.0]]), "NaN"),
        (np.array([[0.0, -1.0], [-1.0, 0.0]]), "negative"),
    ],
)
def test_bad_input_raises_value_error(bad, message):
    with pytest.raises(ValueError, match=message):
        rz.rips(bad)


def test_version_is_exposed():
    assert isinstance(rz.__version__, str)


def test_nan_threshold_raises():
    """A NaN threshold used to admit no edges at all and silently return one
    infinite interval per point."""
    d = rz.pairwise_distances(circle(12))
    with pytest.raises(ValueError, match="nan"):
        rz.rips(d, threshold=math.nan)


def test_degenerate_inputs():
    # a single point: one interval, born at 0, never dies
    (h0,) = rz.rips(np.zeros((1, 1)), maxdim=0)
    assert h0.shape == (1, 2) and math.isinf(h0[0, 1])

    # coincident points merge at distance 0 and are not recorded as bars
    (h0,) = rz.rips(np.zeros((8, 8)), maxdim=0)
    assert h0.shape == (1, 2)

    # asking for more dimensions than there are points is allowed and empty
    d = rz.pairwise_distances(circle(6))
    diagrams = rz.rips(d, maxdim=8)
    assert len(diagrams) == 9
    assert all(x.shape[1] == 2 for x in diagrams)


def test_threshold_below_every_distance_leaves_all_points_separate():
    d = rz.pairwise_distances(circle(10))
    (h0,) = rz.rips(d, maxdim=0, threshold=0.0)
    assert np.isinf(h0[:, 1]).sum() == 10
