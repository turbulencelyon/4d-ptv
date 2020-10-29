import numpy as np

from stm_util import (
    expand_all_neighbours_uniq,
    directional_voxel_traversal,
)

from stm_util_unoptimized import (
    uniquify,
    expand_all_neighbours,
    directional_voxel_traversal2,
)


ps = [(48, 56, 137), (48, 56, 138)]
neighbours = [
    [-1, 0, 0],
    [0, -1, 0],
    [0, 0, -1],
    [0, 0, 1],
    [0, 1, 0],
    [1, 0, 0],
    [0, 0, 0],
]

attended_out_neighbours = [
    [47, 56, 137],
    [48, 55, 137],
    [48, 56, 136],
    [48, 56, 138],
    [48, 57, 137],
    [49, 56, 137],
    [48, 56, 137],
    [47, 56, 138],
    [48, 55, 138],
    [48, 56, 139],
    [48, 57, 138],
    [49, 56, 138],
]


def test_expand_all_neighbours():
    result = uniquify(expand_all_neighbours(ps, neighbours), lambda x: tuple(x))
    assert result == attended_out_neighbours


def test_expand_all_neighbours_uniq():
    result = expand_all_neighbours_uniq(
        np.array(ps, dtype=np.int32), np.array(neighbours, dtype=np.int32)
    )
    tmp = [list(point) for point in result]
    assert sorted(tmp) == sorted(attended_out_neighbours)


point = (-106.10770416259766, -107.7845687866211, 95.57333374023438)
vector_ray = (0.27978962508644883, 0.14909602621096413, 0.9484134861241083)

# fmt: off
cell_bounds = [
    [-140., -126., -112., -98., -84., -70., -56., -42., -28.,
        -14., 0., 14., 28., 42., 56., 70., 84., 98., 112., 126., 140.],
    [-150., -135., -120., -105., -90., -75., -60., -45., -30.,
        -15., 0., 15., 30., 45., 60., 75., 90., 105., 120., 135., 150.],
    [5., 16., 27., 38., 49., 60., 71., 82., 93., 104., 115.,
        126., 137., 148., 159., 170.]
]
attended = [
    (2, 2, 8), (2, 2, 9), (2, 3, 9), (2, 3, 10), (3, 3, 10), (3, 3, 11),
    (3, 3, 12), (3, 3, 13), (3, 3, 14)
]
# fmt: on
cell_bounds_arr = [np.array(bounds) for bounds in cell_bounds]


def test_directional_voxel_traversal2():
    out = directional_voxel_traversal2(point, vector_ray, cell_bounds)
    assert out == attended


def test_directional_voxel_traversal():
    out = directional_voxel_traversal(point, vector_ray, cell_bounds_arr)
    assert np.allclose(out, attended)
