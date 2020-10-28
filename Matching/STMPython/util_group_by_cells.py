# from pprint import pprint

import numpy as np

from transonic import boost
from time import perf_counter


@boost
def special_argsort(cells: "int32[:,:]"):

    max_index = cells.max()
    if max_index > 1024:
        raise ValueError

    multiplicator = 2
    while multiplicator < max_index:
        multiplicator *= 2

    cells = np.array(cells, dtype=np.int32)

    cells_as_ints = (
        multiplicator ** 2 * cells[:, 0]
        + multiplicator * cells[:, 1]
        + cells[:, 2]
    )

    del cells

    indices = cells_as_ints.argsort()

    cells_as_ints_sorted = cells_as_ints[indices]

    return indices, np.array(np.diff(cells_as_ints_sorted), dtype=np.bool)


@boost
def kernel_make_groups_by_cells(
    cam_ray_ids_sorted: "int32[:, :]", diffs: "bool[:]", cam_match: int
):

    groups = []
    start_group = 0
    for index, diff in enumerate(diffs):
        if diff:
            # we exit a group
            stop_group = index + 1
            if stop_group - start_group >= cam_match:
                group = cam_ray_ids_sorted[start_group:stop_group, :]
                cam_ids = group[:, 0]
                nb_cameras = len(set(cam_ids))
                if nb_cameras >= cam_match:
                    groups.append(
                        tuple(
                            tuple(group[index, :])
                            for index in range(group.shape[0])
                        )
                    )
            start_group = stop_group
    return groups


def make_groups_by_cells(cells_all, cam_ray_ids, cam_match: int, log_print=None):
    """

    Inputs
    ------

    traversed:

      Sequence[(cam_id, ray_id, cell)]


    Returns
    -------

    Sequence[Sequence[(cam_id, ray_id)]]

    Set[Tuple[(cam_id, ray_id)]]

    We'd like to filter out groups with less than ``cam_match`` cameras.

    """
    t_start = perf_counter()
    print("PA: make_groups_by_cells")
    log_print(
        "Sorted and grouped by cell index. # of groups:", cells_all.shape[0]
    )

    indices, diffs = special_argsort(cells_all)

    del cells_all

    cam_ray_ids_sorted = cam_ray_ids[indices, :]

    groups = kernel_make_groups_by_cells(cam_ray_ids_sorted, diffs, cam_match)

    print(f"PA # of unique group: {len(groups)}")

    result = tuple(groups)

    print(f"PA: make_groups_by_cells done in {perf_counter() - t_start:.2f} s")

    return result


if __name__ == "__main__":

    cam_ray_ids = np.array([[2, 0], [1, 3], [3, 1], [0, 2]])

    cells = np.array(
        [[1000, 2, 3], [0, 5, 2], [1, 4, 5], [0, 5, 2]], dtype=np.int16
    )

    indices, diffs = special_argsort(cells)

    del cells

    cam_ray_ids_sorted = cam_ray_ids[indices, :]

    groups = kernel_make_groups_by_cells(cam_ray_ids_sorted, diffs, cam_match=2)
