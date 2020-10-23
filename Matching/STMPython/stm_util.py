"""

- 2017-04-12 included [0,0,0] in the default neighbours should not be too
  critical, just a backup measure.

"""

import math
from collections import Counter
import itertools
from copy import copy
import datetime

from typing import Tuple

import numpy as np

from transonic import jit

# import scipy.spatial as sps


def expand_neighbours(p, neighbours):
    "Take a position p and neighbours, create p+n for each n in neighbours"
    # CProfile points towards this line
    return list(list(map(lambda x, y: x + y, p, k)) for k in neighbours)


@jit
def expand_all_neighbours(ps, neighbours):
    """Take positions p and neighbours and create p+n for each n in neighbours
    and for each p

    Note: ~21% of time spent here! (CProfile)
    """
    return joinlists(expand_neighbours(p, neighbours) for p in ps)


def joinlists(lst):
    "'Flatten' a list of lists to a single list"
    # return list(itertools.chain.from_iterable(lst))
    return [value for sub_list in lst for value in sub_list]


def sign(x):
    "Custom sign function that gives either -1 or 1, also for input value 0"
    if x < 0:
        return -1
    else:
        return 1


def special_division(a, b):
    "Division that spawns -infinity in case the denominator is 0"
    if b == 0:
        return -math.inf
    else:
        return a / b


def maplevel(f, item, level):
    "Mathematica style map with level specification"
    if level == 0:
        return f(item)
    else:
        return [maplevel(f, i, level - 1) for i in item]


def uniquify(seq, idfun=None):
    """Return unique elements of a list

    works without second argument for simple elements. Lists-like elements
    require a lambda function lambda x:tuple(x)
    """
    if idfun is None:

        def idfun(x):
            return x

    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen:
            continue
        seen[marker] = 1
        result.append(item)
    return result


def find_index_bin(boundaries, value):
    """
    Gives index in which 'bin' the value value will fall with boundaries 'lst', -1
    = outside

    `boundaries` has to be sorted to make this work properly

    Implemented with O(Log(N)) scaling, basic idea is to reduce by factors of 2
    every time until 1 bin is left

    Greedy: it will take the first bin if the right boundary exactly matches value
    """

    lst = list(boundaries)
    lst[0] -= 10 ** -8
    lst[-1] += 10 ** -8
    mn = 0
    mx = len(lst) - 1
    if lst[mn] <= value <= lst[mx]:
        while mx - mn > 1:
            # print('min',mn,'max',mx)
            # Banker's rounding but does not matter, still O(Log(N)) scaling
            trial = round((mn + mx) / 2)
            if value > lst[trial]:
                mn = trial
            else:
                mx = trial
        # print('Final: min',mn,'max',mx)
        return mn
    else:
        return -1


def vector_norm(y):
    "Euclidean distance of a list (faster than np.linalg.norm(y)!)"
    x = np.array(y)
    return np.sqrt(x.dot(x))


def square_vector_norm(y):
    x = np.array(y)
    return x.dot(x)


def normalize(v):
    norm = vector_norm(v)
    if norm == 0:
        return v
    return v / norm


def closest_point_to_lines2(p, v):
    p1 = np.array(p[0])
    p2 = np.array(p[1])
    v1 = np.array(v[0])
    v2 = np.array(v[1])
    # a = np.dot(v1,v1)          # Assuming v is normalized => length 1
    b = 2 * np.dot(p1 - p2, v1)
    c = 2 * np.dot(v1, v2)
    d = 2 * np.dot(p2 - p1, v2)
    # e = np.dot(v2,v2)          # Assuming v is normalized => length 1
    # f = np.dot(p1, p1) + np.dot(p2, p2)
    # s = (2*a*d + b*c)/(c**2-4*a*e)
    s = (2 * d + b * c) / (c ** 2 - 4)  # Assuming v is normalized => a=e=1
    # t = (c*s - b)/(2*a)
    t = (c * s - b) / 2  # Assuming v is normalized => a=e=1
    sol = (p1 + t * v1 + p2 + s * v2) / 2
    # # Assuming v is normalized => a=1
    # d1 = vector_norm(np.cross(v1,p1-sol))/np.sqrt(a)
    d1 = vector_norm(np.cross(v1, p1 - sol))
    # # Must be the same as d1 for two lines!
    # d2 = np.linalg.norm(np.cross(sol-p1,sol-p1+v1))/np.linalg.norm(v1)
    return [sol.tolist(), d1.item()]


@jit
def step0(a, d):
    length = len(a)
    rhs = np.array([0.0, 0.0, 0.0])
    lhs = length * np.identity(3)
    for i in range(length):
        rhs += a[i] - d[i] * np.dot(a[i], d[i])
        lhs -= np.outer(d[i], d[i])
    return lhs, rhs


# @jit
def step1(a, d, sol):
    dists = list(map(lambda a, d: square_vector_norm(np.cross(sol - a, d)), a, d))
    dists = vector_norm(dists) * np.sqrt(1 / len(a))
    dists = dists.item()  # Convert to regular float
    return [sol.tolist(), dists]


def closest_point_to_lines(p, v):
    """...

    Note: 9% of the time spent here (CProfile)
    """

    if len(p) == 2:
        return closest_point_to_lines2(p, v)
    else:
        a = np.array(p)
        # Assuming v is normalized already
        d = np.array(v)

        # length = len(p)
        # rhs = np.array([0.0, 0.0, 0.0])
        # lhs = length * np.identity(3)
        # for i in range(length):
        #     rhs += a[i] - d[i] * np.dot(a[i], d[i])
        #     lhs -= np.outer(d[i], d[i])

        lhs, rhs = step0(a, d)

        sol = np.linalg.solve(lhs, rhs)

        # dists = list(
        #     map(lambda a, d: square_vector_norm(np.cross(sol - a, d)), a, d)
        # )
        # dists = vector_norm(dists) * np.sqrt(1 / len(a))
        # dists = dists.item()  # Convert to regular float
        # return [sol.tolist(), dists]

        return step1(a, d, sol)


def directional_voxel_traversal2(
    point, vector_direction, cell_bounds, logfile=""
):
    """3D dimensional voxel traversal

    Gives cell-indices back starting from `point` and moving in the direction
    given by `vector_direction`, subject two cell_bounds.

    Note: 12% of time spent here (CProfile)
    """
    assert len(point) == len(vector_direction) == len(cell_bounds) == 3

    # point, vector_direction, and cell_bounds should be 3 dimensional
    if not len(point) == len(vector_direction) == len(cell_bounds) == 3:
        message = "dimension mismatch!"
        if logfile:
            with open(logfile, "a") as flog:
                flog.write(message + "\n")
        raise ValueError(message)
        # print(message)
        # return []

    cell_index_point: Tuple[int, int, int] = tuple(
        map(find_index_bin, cell_bounds, point)
    )

    if -1 in cell_index_point or vector_norm(vector_direction) == 0:
        message = (
            f"ray starts outside cell_bounds!\np: {point}\n"
            f"v: {vector_direction}\ncell index: {cell_index_point}\n"
        )
        if logfile:
            with open(logfile, "a") as flog:
                flog.write(message)

        raise ValueError(message)
        # print(message)
        # return []

    # warning perf: next lines -> could be written with a simple loop (for
    # index_axe in range(3))
    relative_bounds = list(
        map(lambda a, b: [i - b for i in a], cell_bounds, point)
    )
    times = list(
        map(
            lambda a, b: [special_division(i, b) for i in a],
            relative_bounds,
            vector_direction,
        )
    )
    # warning perf: inhomogeneous list [float, int, bool]
    times = list(map(lambda a, b: [[i, b, False] for i in a], times, [0, 1, 2]))
    direction = list(map(sign, vector_direction))
    for index_axe in range(3):
        if direction[index_axe] == 1:
            # vector along this axe (x, y or z)
            times[index_axe][-1][2] = True
        else:
            times[index_axe][0][2] = True
    times = joinlists(times)
    # what is 1 element of times here??
    times = filter(lambda x: x[0] > 0, times)  # Should be > !
    times = sorted(times, key=lambda x: x[0])
    times = list(itertools.takewhile(lambda x: not x[2], times))
    times = [x[1] for x in times]
    out = [copy(cell_index_point)]
    for index in times:
        new = list(out[-1])
        new[index] += direction[index]
        out.append(tuple(new))
    return out


# Dead code... Commented for now
# def directional_voxel_traversal(p, v, bounds, logfile=""):
#     """3D dimensional voxel traversal

#     Gives cell-indices back starting from point p and moving in v direction,
#     subject two cell-bounds bounds
#     """
#     # p, v, and bounds should be 3 dimensional
#     if not len(p) == len(v) == len(bounds) == 3:
#         print("dimension mismatch!")
#         if logfile:
#             with open(logfile, "a") as flog:
#                 flog.write("dimension mismatch!\n")
#         return []

#     curpos = list(p)
#     cellindex = list(map(find_index_bin, bounds, p))
#     if -1 in cellindex or vector_norm(v) == 0:
#         message = (
#             "ray starts outside bounds!\np: {p}\nv: {v}\n"
#             "cell index: {cellindex}\n"
#         )
#         if logfile:
#             with open(logfile, "a") as flog:
#                 flog.write(message)

#         print(message)
#         return []

#     sgns = list(map(sign, v))
#     sgnspart = list(map(lambda x: 1 if x == -1 else 0, sgns))

#     cont = True
#     steps = 0
#     out = [copy(cellindex)]
#     while cont and steps < 10000:
#         # Steps is just a safety thing for now, can be removed later
#         steps += 1
#         newindices = list(map(lambda x, y: x + y, cellindex, sgns))
#         newbounds = [0, 0, 0]
#         for i in range(3):
#             if 0 <= newindices[i] <= len(bounds[i]) - 1:
#                 # Should there be + sgnspart[i] in the middle term?
#                 newbounds[i] = bounds[i][newindices[i] + sgnspart[i]]
#             else:
#                 # Bounds are at infinity, so the 'next' bounds are infinitely far away
#                 newbounds[i] = math.inf

#         if math.inf not in newbounds:
#             ts = [0, 0, 0]
#             for i in range(3):
#                 if v[i] == 0:
#                     ts[i] = math.inf  # It will take infinite amount of time
#                 else:
#                     ts[i] = (newbounds[i] - curpos[i]) / v[i]

#             # Find the 'times' needed to the next boundaries in each dimensions
#             order = sorted(range(len(ts)), key=lambda k: ts[k])
#             # Find the dimension for which the time is shortest
#             minpos = order[0]

#             cellindex[minpos] += sgns[minpos]
#             ts = ts[minpos]
#             for i in range(3):
#                 curpos[i] += v[i] * ts

#             out.append(copy(cellindex))
#         else:
#             # print("at edge")
#             cont = False
#     return out


def at_face(bmin, bmax, hitb):
    return bmin <= hitb <= bmax


def prepare_ray(p, v, bounds):
    """
    Projects ray (defined by p, v) onto an AABB (axis aligned bounding box).

    Returns (bool_hit, bool_inside, position, vector_direction)

    - bool_hit tells if it hits AABB,
    - bool_inside tells if it is projected
    - position: the new position or [] in case it misses.
    - vector_direction: the normalized vector

    Note that ray can be projected onto an AABB with negative 'time'...

    Note: very quick so no need to optimize.
    """

    xmin = bounds[0][0]
    xmax = bounds[0][1]
    ymin = bounds[1][0]
    ymax = bounds[1][1]
    zmin = bounds[2][0]
    zmax = bounds[2][1]
    x = p[0]
    y = p[1]
    z = p[2]
    vector_direction = normalize(v).tolist()  # v is normalized
    vx, vy, vz = vector_direction
    if xmin < x < xmax and ymin < y < ymax and zmin < z < zmax:
        return [
            True,
            True,
            p,
            vector_direction,
        ]  # Return False and original point
    else:
        t = list(
            map(
                special_division,
                [xmin - x, xmax - x, ymin - y, ymax - y, zmin - z, zmax - z],
                [vx, vx, vy, vy, vz, vz],
            )
        )
        ip = [0 for tmp in range(6)]
        for i in range(6):
            ti = t[i]
            if abs(ti) == math.inf:
                ip[i] = [math.inf, math.inf, math.inf]
            else:
                ip[i] = [x + vx * ti, y + vy * ti, z + vz * ti]

    atfacex1 = at_face(ymin, ymax, ip[0][1]) and at_face(zmin, zmax, ip[0][2])
    atfacex2 = at_face(ymin, ymax, ip[1][1]) and at_face(zmin, zmax, ip[1][2])
    atfacey1 = at_face(xmin, xmax, ip[2][0]) and at_face(zmin, zmax, ip[2][2])
    atfacey2 = at_face(xmin, xmax, ip[3][0]) and at_face(zmin, zmax, ip[3][2])
    atfacez1 = at_face(xmin, xmax, ip[4][0]) and at_face(ymin, ymax, ip[4][1])
    atfacez2 = at_face(xmin, xmax, ip[5][0]) and at_face(ymin, ymax, ip[5][1])
    data = [t, [atfacex1, atfacex2, atfacey1, atfacey2, atfacez1, atfacez2], ip]
    # Data will be a list of lists, each having the form: time-till-hit, hits
    # face (boolean), position it hits a plane
    data = list(zip(*data))
    # Select only those that hit a face #don't change to 'is True' numpy bool
    # possibility
    data = list(filter(lambda xx: xx[1], data))
    if len(data) > 0:
        # Sort by arrival time (time till hit)
        data = sorted(data, key=lambda x: x[0])
        # Position it hits the plane of first-hit
        return [True, False, data[0][2], vector_direction]
    else:
        return [False, False, [], vector_direction]


def space_traversal_matching(
    raydata,
    boundingbox,
    nx=75,
    ny=75,
    nz=75,
    cam_match_func=lambda x: len(x) > 2,
    max_matches_per_ray=2,
    maxdistance=999.9,
    neighbours=6,
    logfile="",
):

    if not (
        len(boundingbox) == 3
        and len(boundingbox[0]) == 2
        and len(boundingbox[1]) == 2
        and len(boundingbox[2]) == 2
        and nx >= 5
        and ny >= 5
        and nz >= 5
    ):
        print(
            "Something went wrong:\n"
            "Bounding box should be of the form "
            "[[xmin,xmax],[ymin,ymax],[zmin,zmax]] {boundingbox}\n",
            "nx,ny,nz should be >=5 otherwise you will do a lot of comparisons!",
            [nx, ny, nz],
        )
        return []

    if neighbours == 0:
        neighbours = [[0, 0, 0]]
    elif neighbours == 6:
        neighbours = [
            [-1, 0, 0],
            [0, -1, 0],
            [0, 0, -1],
            [0, 0, 1],
            [0, 1, 0],
            [1, 0, 0],
            [0, 0, 0],
        ]
    elif neighbours == 18:
        neighbours = [
            [-1, -1, 0],
            [-1, 0, -1],
            [-1, 0, 0],
            [-1, 0, 1],
            [-1, 1, 0],
            [0, -1, -1],
            [0, -1, 0],
            [0, -1, 1],
            [0, 0, -1],
            [0, 0, 1],
            [0, 1, -1],
            [0, 1, 0],
            [0, 1, 1],
            [1, -1, 0],
            [1, 0, -1],
            [1, 0, 0],
            [1, 0, 1],
            [1, 1, 0],
            [0, 0, 0],
        ]
    elif neighbours == 26:
        neighbours = [
            [-1, -1, -1],
            [-1, -1, 0],
            [-1, -1, 1],
            [-1, 0, -1],
            [-1, 0, 0],
            [-1, 0, 1],
            [-1, 1, -1],
            [-1, 1, 0],
            [-1, 1, 1],
            [0, -1, -1],
            [0, -1, 0],
            [0, -1, 1],
            [0, 0, -1],
            [0, 0, 1],
            [0, 1, -1],
            [0, 1, 0],
            [0, 1, 1],
            [1, -1, -1],
            [1, -1, 0],
            [1, -1, 1],
            [1, 0, -1],
            [1, 0, 0],
            [1, 0, 1],
            [1, 1, -1],
            [1, 1, 0],
            [1, 1, 1],
            [0, 0, 0],
        ]
    if not isinstance(neighbours, list):
        raise TypeError(
            "Neighbours should be a 0, 6, 18, or 26 or "
            f"a list of triplets: {neighbours}"
        )

    if logfile:

        def log_print(*args):
            with open(logfile, "a") as flog:
                flog.write(" ".join(map(str, list(args))) + "\n")

    else:

        def log_print(*args):
            print(*args)

    log_print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    bounds = list(
        map(
            lambda x, n: [x[0] + i * (x[1] - x[0]) / n for i in range(n + 1)],
            boundingbox,
            [nx, ny, nz],
        )
    )

    bounds_pa = [
        list(np.linspace(limits[0], limits[1], n + 1))
        for limits, n in zip(boundingbox, (nx, ny, nz))
    ]

    assert np.allclose(bounds[0], bounds_pa[0]), (
        np.array(bounds_pa[0]).dtype,
        np.array(bounds[0]).dtype,
    )

    for i in range(3):
        assert bounds[i][-1] == boundingbox[i][-1]
        assert bounds_pa[i][-1] == boundingbox[i][-1]

    log_print("# of cells:", nx * ny * nz)
    rays = list(raydata)
    print(len(raydata), raydata[0])

    # Prepare the rays so that they are inside the box or on the side.
    # First element is camera ID
    INDEX_CAM = 0
    # Second element is the ray ID
    INDEX_RAY = 1

    def cam_marker_func(x):
        return x[INDEX_CAM]

    def ray_marker_func(x):
        return x[INDEX_RAY]

    # Store a dictionary of (cameraid, ray_id): [pos, direction]
    raydb = {}
    valid_rays = []  # Store the transformed rays
    numrays = Counter()  # Store the number of rays per camera
    invalidcounter = Counter()  # Store the number that are invalid
    for ray in rays:
        cam_id = ray[INDEX_CAM]
        ray_id = ray[INDEX_RAY]
        numrays[cam_id] += 1

        pp = list(ray[2:5])
        vv = list(ray[5:8])
        bool_hit, bool_inside, position, vector_direction = prepare_ray(
            pp, vv, boundingbox
        )

        if bool_hit:  # If it does not miss (hit or inside)
            raydb[(cam_id, ray_id)] = [position, vector_direction]
            valid_rays.append(
                [cam_id, ray_id, bool_inside, position, vector_direction]
            )
        else:
            invalidcounter[cam_id] += 1

    # for k, v in raydb.items():
    #    print(vector_norm(v[1]))

    log_print(
        "# of rays for each camera:",
        dict(numrays),
        "\n# of rays that miss the bounding box:",
        dict(invalidcounter),
    )
    if len(dict(numrays)) > 10:
        log_print(
            "# of cameras is large:",
            len(dict(numrays)),
            " Double check the input!",
        )
        return []

    # traversed: List[List[int, int, List[int, int, int]]]
    # [cam_id, ray_id, cell]
    traversed = []
    for ray in valid_rays:
        cam_id, ray_id, bool_inside, position, vector_direction = ray
        if bool_inside:  # Ray is inside, traverse both forward and backward
            out = directional_voxel_traversal2(
                position, vector_direction, bounds, logfile
            ) + directional_voxel_traversal2(
                position,
                list(map(lambda x: -1 * x, vector_direction)),
                bounds,
                logfile,
            )
        else:  # Ray is at edge, traverse in forward direction only
            out = directional_voxel_traversal2(
                position, vector_direction, bounds, logfile
            )

        # Expand in neighbourhood and remove duplicates
        cells = uniquify(
            expand_all_neighbours(out, neighbours), lambda x: tuple(x)
        )
        # Creates long list of cam_id, ray_id, cellindex
        ext = [[cam_id, ray_id, cell] for cell in cells]
        traversed.extend(ext)  # Pile up these into traversed

    log_print("# of voxels traversed after expansion:", len(traversed))

    def cell_func(x):
        return x[2]

    # warning perf: this could be rewritten and accelerated with Pythran?
    # Sort based on cell. Sort is needed before groupby
    # (CProfile points towards this line)
    traversed = sorted(traversed, key=cell_func)
    # Group elements by same cell (CProfile points towards this line)
    traversed = [list(g) for k, g in itertools.groupby(traversed, cell_func)]
    log_print("Sorted and grouped by cell index. # of groups:", len(traversed))
    # Prune based on number of rays (fast rough filter, cam filter later)
    traversed = list(filter(cam_match_func, traversed))
    log_print("Rough pruned based on number of cameras:", len(traversed))
    # Remove cellindex, not needed anymore, leave [cam_id, ray_id]
    traversed = maplevel(lambda x: [x[0], x[1]], traversed, 2)
    traversed = list(
        map(
            lambda x: [list(g) for k, g in itertools.groupby(x, cam_marker_func)],
            traversed,
        )
    )
    # log_print("Cell index removed and grouped by camera for each cell")
    # Prune based on number of different cameras
    traversed = tuple(filter(cam_match_func, traversed))
    log_print("Pruned based on number of cameras:", len(traversed))
    # All combinations between all cameras
    candidates = list(
        map(
            lambda x: [list(tup) for tup in itertools.product(*x)],
            traversed,
        )
    )
    # Flatten a list of lists to a single list
    candidates = joinlists(candidates)
    log_print("Flattened list of candidates:", len(candidates))
    # Delete duplicates, flattened list as tag
    candidates = uniquify(candidates, lambda x: tuple(joinlists(x)))
    log_print("Duplicate candidates removed:", len(candidates))
    candidates = sorted(candidates)

    # log_print("Computing match position and quality of candidates...")
    newcandidates = []
    for c in candidates:
        pvdata = [raydb[tuple(x)] for x in c]
        pdata = [x[0] for x in pvdata]
        vdata = [x[1] for x in pvdata]
        out = closest_point_to_lines(pdata, vdata)
        newcandidates.append([c, list(out[0]), out[1]])

    # hullpts = [[20.0,5.0,165.0],[20.0,10.0,160.0],[20.0,10.0,165.0],[20.0,15.0,160.0],[20.0,15.0,165.0],[20.0,20.0,160.0],[20.0,20.0,165.0],[20.0,25.0,160.0],[20.0,25.0,165.0],[25.0,0.0,165.0],[25.0,0.0,170.0],[25.0,5.0,155.0],[25.0,10.0,155.0],[25.0,15.0,155.0],[25.0,20.0,155.0],[25.0,25.0,155.0],[25.0,25.0,175.0],[25.0,30.0,170.0],[25.0,35.0,155.0],[25.0,35.0,160.0],[30.0,-5.0,170.0],[30.0,0.0,160.0],[30.0,5.0,150.0],[30.0,10.0,150.0],[30.0,15.0,150.0],[30.0,20.0,150.0],[30.0,25.0,150.0],[30.0,25.0,180.0],[30.0,30.0,150.0],[30.0,30.0,180.0],[30.0,35.0,175.0],[30.0,40.0,165.0],[30.0,45.0,155.0],[35.0,-10.0,175.0],[35.0,-5.0,165.0],[35.0,-5.0,180.0],[35.0,0.0,155.0],[35.0,10.0,145.0],[35.0,15.0,145.0],[35.0,20.0,145.0],[35.0,25.0,185.0],[35.0,30.0,185.0],[35.0,35.0,150.0],[35.0,35.0,185.0],[35.0,40.0,180.0],[35.0,45.0,155.0],[35.0,45.0,170.0],[35.0,50.0,160.0],[40.0,-10.0,175.0],[40.0,-10.0,180.0],[40.0,-5.0,165.0],[40.0,0.0,155.0],[40.0,5.0,145.0],[40.0,10.0,145.0],[40.0,15.0,145.0],[40.0,30.0,150.0],[40.0,40.0,180.0],[40.0,45.0,155.0],[40.0,45.0,175.0],[40.0,50.0,160.0],[40.0,50.0,165.0],[45.0,-5.0,175.0],[45.0,0.0,165.0],[45.0,5.0,155.0],[45.0,10.0,150.0],[45.0,15.0,150.0],[45.0,40.0,180.0],[45.0,50.0,160.0],[45.0,50.0,165.0],[50.0,15.0,160.0],[50.0,20.0,160.0],[50.0,25.0,175.0],[50.0,30.0,175.0],[50.0,35.0,175.0],[50.0,45.0,165.0],[55.0,15.0,170.0],[55.0,20.0,170.0],[55.0,25.0,170.0],[55.0,30.0,170.0],[55.0,35.0,170.0],[55.0,40.0,170.0]];
    # delaun = sps.Delaunay(hullpts)  # Define the Delaunay triangulation
    # inq = delaun.find_simplex([x[1] for x in newcandidates])>0
    # inq = inq.tolist();
    # print(inq)
    # print(type(inq))
    # newcandidates = list(map(lambda x,y: x + [y],newcandidates,inq))

    candidates = copy(newcandidates)
    del newcandidates
    # log_print("Sorting candidate matches by quality of match...")
    # sort by number of cameras then by error
    candidates = sorted(candidates, key=lambda x: (-len(x[0]), x[2]))
    # log_print("Here are upto 9999 of the best matches:")
    # log_print("Index, [cam_id ray_id ....] Position, Mean square distance")
    # print("num candidates:",len(candidates))
    # for i in range(1,len(candidates),100):
    #    print(i,candidates[i])

    log_print(
        f"Selecting the best matches with up to {max_matches_per_ray}",
        f"match(es)/ray out of {len(candidates)} candidates",
    )

    # now we want to pick the best matches first and match each ray at most
    # max_matches_per_ray

    approvedmatches = []  # Store approved candidates
    matchcounter = Counter()  # Keep track of how many they are matched
    for cand in candidates:
        if cand[2] < maxdistance:
            valid = True
            for idpair in cand[0]:
                if matchcounter[tuple(idpair)] >= max_matches_per_ray:
                    valid = False
                    # print(idpair,"has been matched already",max_matches_per_ray,"time(s)")
                    break
            if valid:
                for idpair in cand[0]:
                    matchcounter[tuple(idpair)] += 1

                approvedmatches.append(list(cand))

    log_print(
        "Selecting done.",
        len(approvedmatches),
        "matched found (out of",
        len(candidates),
        "candidates)",
    )
    # print("Here are the approved matches:")
    # print("Index, [cam_id ray_id ....] Position, Mean square distance")

    return approvedmatches
