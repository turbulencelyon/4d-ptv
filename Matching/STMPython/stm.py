#! /usr/bin/env python
"""
Space Traversal Matching: compute matches from rays projecting them into voxels.

Example:

  export PATH_INPUT_DATA="../../Documentation/TestData/Processed_DATA/MyExperiment/Parallel/Matching/Rays/rays_1-10.dat"
  ./stm.py $PATH_INPUT_DATA 1 2 2 0.2 400 400 250 2

or (to specify the limits of the visualized region):

  ./stm.py $PATH_INPUT_DATA 1 2 2 0.2 400 400 250 2 "[[-140, 140], [-150, 150], [5, 170]]"

"""
import os
import copy
import struct
import traceback
from time import perf_counter
import argparse


USE_UNOPTIMIZED = os.environ.get("STM_PYTHON_USE_UNOPTIMIZED", False)

"""
For xonsh::

  $PATH_INPUT_DATA="../../Documentation/TestData/Processed_DATA/MyExperiment/Parallel/Matching/Rays/rays_1-10.dat"

Note: the "unoptimized" code is actually slightly faster than the first
Python code: 179.5 s versus 184 s, i.e. 2.4 % faster :-)

To be compared to the perf of the C++ code: 8.45 s!

- The old Python code was 21.8 time slower than the optimized C++ one.

- The unoptimized Python code is 21.2 time slower than the optimized C++ one.

Command to launch the same computation (?) with Python::

  STM_PYTHON_USE_UNOPTIMIZED=1 python stm.py $PATH_INPUT_DATA 1 2 2 0.2 400 400 250 2
  python stm.py $PATH_INPUT_DATA 1 2 2 0.2 400 400 250 2

With both versions, the same number of matches are found out of the same number
of candidates: "6754 matched found (out of 135271 candidates)".

Command to launch the C++ code::

  ../STMCpp/STM -i $PATH_INPUT_DATA -f 1 -c 2 -d 0.2 -m 2 -x 400 -y 400 -z 250 -b -140 140 -150 150 5 170 --hdf5

The result of the C++ code is 7466 matched found (out of 135271 candidates)
"""

if not USE_UNOPTIMIZED:
    from stm_util import space_traversal_matching
else:
    print("Using stm_util_unoptimized")
    from stm_util_unoptimized import space_traversal_matching


def compute_stm(
    filename,
    start_frame,
    stop_frames,
    cam_match,
    max_distance,
    nx,
    ny,
    nz,
    max_matches_per_ray,
    bounding_box=[[-140, 140], [-150, 150], [5, 170]],
    min_distance_matches_1ray=None,
    neighbours=6,
):
    """
    Compute matches from rays projecting them into voxels.

    See the output of ``./stm.py -h`` for the meaning of the arguments.

    - neighbours

    Number of illuminated voxels: due to noise, when a ray crosses a voxel, it
    is possible that in reality, the ray crosses a close voxel. neighbours
    indicates how many neighbours we consider in reality when a ray crosses a
    voxel. =6 by defaut.

    """
    #############################################################################################################
    # Parameters to adjust
    tstart = perf_counter()

    #############################################################################################################

    fileout = copy.copy(filename).split(".")
    fileout = ".".join(fileout[0 : len(fileout) - 1])
    filelog = fileout + ".log"
    fileout = (
        fileout.replace("rays", "matched")
        + f"cam{cam_match}_{start_frame}-{stop_frames-1}.dat"
    )

    fout = open(fileout, "wb")
    fin = open(filename, "rb")
    frameid = start_frame
    numpts = fin.read(4)  # Read 4 bytes header
    while len(numpts) > 0 and frameid < stop_frames:  # If something is read
        numpts = struct.unpack("I", numpts)[0]  # Interpret header as 4 byte uint
        with open(filelog, "a") as flog:
            flog.write("#######\nFrame: {frameid}\nNumber of rays: {numpts}\n")

        print("Frame:", frameid, ". # of rays:", numpts)

        # Read rays
        raydata = fin.read(numpts * 27)  # 27 bytes per line 2+1+6*4
        raydata = struct.unpack(
            "=" + ("BH6f" * numpts), raydata
        )  # Create string '=BHFFFFFFBHFFFFFFBHFFFFFF...BHFFFFFF'
        raydata = list(
            map(
                lambda i: list(raydata[8 * i : 8 * i + 8]),
                range(len(raydata) // 8),
            )
        )  # Reshape to 8*N np.arreyreshape converts everything to floats...
        # The actual call
        try:
            output = space_traversal_matching(
                raydata,
                bounding_box,
                nx=nx,
                nz=nz,
                ny=ny,
                cam_match=cam_match,
                neighbours=neighbours,
                logfile=filelog,
                max_distance=max_distance,
                min_distance_matches_1ray=min_distance_matches_1ray,
            )
        except ValueError:
            tb = traceback.format_exc()
            with open(filelog, "a") as flog:
                flog.write(tb + "\n")

            raise

        # Prepare output
        print("Matches found:", len(output))
        coutput = []
        for o in output:
            tmp = [len(o[0])]
            tmp.extend(o[1])
            tmp.append(o[2])
            for j in o[0]:
                tmp.extend(j)
            coutput.append(tmp)
        output = coutput
        del coutput

        # Output the output
        buf = struct.pack("I", len(output))  # Write the number of matches
        fout.write(buf)

        for out in output:
            buf = struct.pack(
                "=B4f" + out[0] * "BH", *out
            )  # Write each to the output file
            fout.write(buf)

        numpts = fin.read(4)  # Read next header
        frameid += 1
    fout.close()
    fin.close()
    print("Finished")

    elapsed = perf_counter() - tstart
    print(f"Elapsed time:       {elapsed:.2f} s")
    print(f"Elapsed time/frame: {elapsed / (stop_frames - start_frame):.2f} s")


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "path_file",
        type=str,
        help="Path towards the file containing the ray data",
    )

    parser.add_argument(
        "start_frame",
        type=int,
        help="Index of the first frame",
    )

    parser.add_argument(
        "stop_frame",
        type=int,
        help="Index of the last frame + 1",
    )

    parser.add_argument(
        "cam_match",
        type=int,
        help="Minimum number of rays crossing to get a match",
    )

    parser.add_argument(
        "max_distance",
        type=float,
        help="Maximum distance allowed for a match",
    )

    parser.add_argument(
        "nx",
        type=int,
        help="Number of voxels in the x direction",
    )

    parser.add_argument(
        "ny",
        type=int,
        help="Number of voxels in the y direction",
    )

    parser.add_argument(
        "nz",
        type=int,
        help="Number of voxels in the z direction",
    )

    parser.add_argument(
        "max_matches_per_ray",
        type=int,
        help="Maximum number of matches/ray",
    )

    parser.add_argument(
        "bounding_box",
        type=str,
        nargs="?",
        help=(
            "Corresponds to the volume visualized "
            "[[minX, maxX], [minY, maxY], [minZ, maxZ]]"
        ),
        default="[[-140, 140], [-150, 150], [5, 170]]",
    )

    parser.add_argument(
        "-md1r",
        "--min-distance-matches-1ray",
        type=float,
        default=None,
        help="Minimum distance for multiple matches per ray",
    )

    args = parser.parse_args()
    bounding_box_as_str = args.bounding_box
    args.bounding_box = eval(bounding_box_as_str)

    return args


def main():

    args = parse_args()
    print(args)

    compute_stm(
        args.path_file,
        args.start_frame,
        args.stop_frame,
        args.cam_match,
        args.max_distance,
        args.nx,
        args.ny,
        args.nz,
        args.max_matches_per_ray,
        args.bounding_box,
        args.min_distance_matches_1ray,
    )


if __name__ == "__main__":
    main()
