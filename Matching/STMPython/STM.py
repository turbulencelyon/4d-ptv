#! /usr/bin/env python
"""

Example::

  python STM.py "../../Documentation/TestData/Processed_DATA/MyExperiment/Parallel/Matching/Rays/rays_1-10.dat" 1 2 2 0.2 400 400 250 2

"""
import socket
import sys
import copy
import struct
from datetime import datetime

import numpy as np

from STMFunctions import space_traversal_matching

print(socket.gethostname())
print("Python", sys.version)


def STM(
    filename,
    minframes,
    maxframes,
    cam_match,
    maxdistance,
    nx,
    ny,
    nz,
    max_matches_per_ray,
    boundingbox=[[-140, 140], [-150, 150], [5, 170]],
    neighbours=6,
):
    """
    Compute matches from rays projecting them into voxels.

    # Parameters
    #   filename                           # name of the file containing rays
    #   minframes                          # number of the first frame
    #   maxframes                          # number of the last frame
    #   cam_match                           # minimum number of rays crossing to get a match. We require a match to satisfy cam_match_func = lambda x: len(x)>=cam_match for the number of cameras (and thus the number of rays)
    #   maxdistance                        # max distance allowed for a match.
    #   nx,ny,nz                           # number of voxels in each direction
    #   max_matches_per_ray                   # number of matches/ray
    #   boundingbox                        # correspond to the volume visualized [[minX,maxX],[minY,maxY],[minZ,maxZ]] ATTENTION Does not work currently -> To DO !!
    #   neighbours                         # number of illuminated voxels: due to noise, when a ray crosses a voxel, it is possible that in reality, the ray crosses a close voxel. neighbours indicates how many neighbours we consider in reality when a ray crosses a voxel. =6 by defaut.
    """
    #############################################################################################################
    # Parameters to adjust
    tstart = datetime.now().timestamp()
    cam_match_func = (
        lambda x: len(x) >= cam_match
    )  # We require a match to satisfy this requirement for the number of cameras (and thus the number of rays)
    #############################################################################################################

    fileout = copy.copy(filename).split(".")
    fileout = ".".join(fileout[0 : len(fileout) - 1])
    filelog = fileout + ".log"
    print(filelog)
    fileout = (
        fileout.replace("rays", "matched")
        + f"cam{cam_match}_{minframes}-{maxframes}.dat"
    )

    fout = open(fileout, "wb")
    fin = open(filename, "rb")
    frameid = minframes
    numpts = fin.read(4)  # Read 4 bytes header
    while len(numpts) > 0 and frameid < maxframes:  # If something is read
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
        output = space_traversal_matching(
            list(raydata),
            boundingbox,
            nx=nx,
            nz=nz,
            ny=ny,
            cam_match_func=cam_match_func,
            neighbours=neighbours,
            logfile=filelog,
            maxdistance=maxdistance,
        )

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

    elapsed = datetime.now().timestamp() - tstart
    print("Elapsed time:", elapsed)
    print("Elapsed time/frame:", elapsed / (frameid - 1))


if len(sys.argv) == 10:
    STM(
        str(sys.argv[1]),
        int(sys.argv[2]),
        int(sys.argv[3]),
        int(sys.argv[4]),
        float(sys.argv[5]),
        int(sys.argv[6]),
        int(sys.argv[7]),
        int(sys.argv[8]),
        int(sys.argv[9]),
    )
elif len(sys.argv) == 11:
    STM(
        str(sys.argv[1]),
        int(sys.argv[2]),
        int(sys.argv[3]),
        int(sys.argv[4]),
        float(sys.argv[5]),
        int(sys.argv[6]),
        int(sys.argv[7]),
        int(sys.argv[8]),
        int(sys.argv[9]),
        np.array(sys.argv[10]),
    )
else:
    print(
        f"Only {len(sys.argv) - 1} arguments\n"
        "There should be an argument with the "
        "filename, minframe, maxframe, cam_match, maxdistance, nx, ny, nz, "
        "max_matches_per_ray, boundingbox(optional)!"
    )
