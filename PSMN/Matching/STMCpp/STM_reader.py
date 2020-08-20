"""STM_reader
=============

This module implements reading the original STM binary files, and the STM
HDF5 files.
"""

from typing import List
from dataclasses import dataclass
import sys
import h5py
import numpy as np


def fromfile(fid, dtype, count):
    return np.frombuffer(fid.read(np.dtype(dtype).itemsize * count), dtype, count)


@dataclass
class Camrayid:
    """This class represents a pair <camid, rayid>.
    It is similar to the C++ camrayid structure.
    """

    camid: int
    rayid: int


@dataclass
class CandidateMatch:
    """This class represents a match.
    It is similar to the C++ candidatematch structure.
    """

    x: float
    y: float
    z: float
    error: float
    camrayid: List[Camrayid]


class STM_Reader:
    """This class implements the reading of STM files.
    """

    def __init__(self, filename):

        self.filename = filename
        self.frames = list()
        if filename.endswith(".h5") or filename.endswith(".hdf5"):
            self.fmt = "hdf5"
            self.read_hdf5()
        elif filename.endswith(".bin"):
            self.fmt = "bin"
            self.read_bin()
        else:
            raise RuntimeError("Unknown file format")

    def read_hdf5(self, verbose=True):
        with h5py.File(self.filename, "r") as f:
            self.maxframes, = f.attrs["maxframes"]
            self.mincameras, = f.attrs["mincameras"]
            self.maxdistance, = f.attrs["maxdistance"]
            self.maxmatchesperray, = f.attrs["maxmatchesperray"]
            self.nx, = f.attrs["nx"]
            self.ny, = f.attrs["ny"]
            self.nz, = f.attrs["nz"]

            if verbose:
                print("maxframes =", self.maxframes)
                print("mincameras =", self.mincameras)
                print("maxdistance =", self.maxdistance)
                print("maxmatchesperray =", self.maxmatchesperray)
                print("nx =", self.nx)
                print("ny =", self.ny)
                print("nz =", self.nz)

            for frameno in range(1, self.maxframes + 1):
                xyze_data = f[f"frame{frameno:}_xyze"][()]
                camrayids_data = f[f"frame{frameno:}_camrayids"][()]
                if verbose:
                    print(f"Reading Frame {frameno:} datasets")
                n, nrows = xyze_data.shape
                if n != 4:
                    raise TypeError("Wrong size for xyze dataset")
                m, nrows_bis = camrayids_data.shape
                if nrows != nrows_bis:
                    raise TypeError("Size mismatch between xyze and camrayids datasets")

                matches = list()
                for row in range(nrows):
                    camids = camrayids_data[0::2, row]
                    rayids = camrayids_data[1::2, row]
                    camrayids = [
                        Camrayid(camid, rayid)
                        for camid, rayid in zip(camids, rayids)
                        if camid != -1 and rayid != -1
                    ]
                    matches.append(
                        CandidateMatch(
                            xyze_data[0, row],
                            xyze_data[1, row],
                            xyze_data[2, row],
                            xyze_data[3, row],
                            camrayids,
                        )
                    )
                self.frames.append(matches)

    def read_bin(self):
        with open(self.filename, 'rb') as f:
            currentframe = 0
            while True:
                try:
                    numberofmatches, = fromfile(f, "<u4", 1)
                except ValueError:
                    break
                if numberofmatches:
                    currentframe += 1
                    print(f"Reading data for frame {currentframe:}")
                    matches = list()
                    for matchno in range(numberofmatches):
                        numberofcams, = fromfile(f, "B", 1)
                        x, y, z, err = fromfile(f, "<f4", 4)
                        camrayids = list()
                        for idno in range(numberofcams):
                            camid, = fromfile(f, "B", 1)
                            rayid, = fromfile(f, "<u2", 1)
                            camrayids.append(Camrayid(camid, rayid))
                        matches.append(CandidateMatch(x, y, z, err, camrayids))
                    self.frames.append(matches)
                else:
                    break

    def describe(self):
        for frameno, frame in enumerate(self.frames):
            print(f"Frame {frameno:}")
            print("    No matches:", len(frame))


if __name__ == '__main__':
    basename = sys.argv[1]
    if basename.endswith("_out_cpp.h5"):
        basename = basename[:-11]

    bin_file = STM_Reader(basename + "_out_cpp.bin")
    print("Bin file")
    print("========")
    bin_file.describe()
    print()

    h5_file = STM_Reader(basename + "_out_cpp.h5")
    print("HDF5 file")
    print("=========")
    h5_file.describe()
    print()

    print("Comparing files")
    similar = True
    for frame_h5, frame_bin in zip(h5_file.frames, bin_file.frames):
        for match_h5, match_bin in zip(frame_h5, frame_bin):
            # à cause de la différence entre float et double, on ne peut pas tester l'égalité
            # if match_h5 != match_bin:
            #     print("Difference found !")
            #     print("match_h5 =", match_h5)
            #     print("match_bin =", match_bin)
            #     similar = False
            #     break
            h5_data = np.array([match_h5.x, match_h5.y, match_h5.z, match_h5.error])
            bin_data = np.array([match_bin.x, match_bin.y, match_bin.z, match_bin.error])
            if not np.isclose(h5_data, bin_data).all():
                print("Difference found in (x, y, z, err) !")
                print("match_h5 =", match_h5)
                print("match_bin =", match_bin)
                similar = False
                break
            if match_h5.camrayid != match_bin.camrayid:
                print("Difference found in camrayid !")
                print("match_h5 =", match_h5)
                print("match_bin =", match_bin)
                similar = False
                break

        if not similar:
            break

    if similar:
        print("Files are identical")
    else:
        print("Files are different")
