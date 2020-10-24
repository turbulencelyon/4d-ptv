from stm_util import expand_all_neighbours


def test_expand_all_neighbours():
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

    result = expand_all_neighbours(ps, neighbours)

    attended = [
        [47, 56, 137],
        [48, 55, 137],
        [48, 56, 136],
        [48, 56, 138],
        [48, 57, 137],
        [49, 56, 137],
        [48, 56, 137],
        [47, 56, 138],
        [48, 55, 138],
        [48, 56, 137],
        [48, 56, 139],
        [48, 57, 138],
        [49, 56, 138],
        [48, 56, 138],
    ]

    assert result == attended
