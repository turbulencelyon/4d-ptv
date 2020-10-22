"""
Example of output:

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1   40.693   40.693   44.036   44.036 STMPython/stm_util.py:576(<listcomp>)
 26944128   29.400    0.000   36.194    0.000 STMPython/stm_util.py:21(<genexpr>)
   549354   18.742    0.000   21.074    0.000 {built-in method builtins.sorted}
    61386   15.341    0.000   15.341    0.000 STMPython/stm_util.py:216(<listcomp>)
        1   13.056   13.056  232.238  232.238 STMPython/stm_util.py:377(space_traversal_matching)
    10429   10.052    0.001   17.929    0.002 STMPython/stm_util.py:59(uniquify)
  2568545    7.653    0.000    8.168    0.000 STMPython/stm_util.py:585(<listcomp>)
    20462    7.398    0.000   38.002    0.002 STMPython/stm_util.py:176(directional_voxel_traversal2)
 70728336    6.794    0.000    6.794    0.000 STMPython/stm_util.py:21(<lambda>)
  1056990    5.115    0.000    6.140    0.000 .../python3.8/site-packages/numpy/core/numeric.py:1308(normalize_axis_tuple)
    10428    4.864    0.000    4.864    0.000 STMPython/stm_util.py:568(<listcomp>)
  3368016    4.809    0.000   41.003    0.000 STMPython/stm_util.py:19(expand_neighbours)

"""

import pstats
import cProfile
from time import perf_counter

from stm import compute_stm

args = [
    "../../Documentation/TestData/Processed_DATA/"
    "MyExperiment/Parallel/Matching/Rays/rays_1-10.dat"
]

args.extend([1, 2, 2, 0.2, 400, 400, 250, 2])


t0 = perf_counter()
cProfile.runctx("compute_stm(*args)", globals(), locals(), "profile.pstats")

s = pstats.Stats("profile.pstats")
s.sort_stats("time").print_stats(12)

t_end = perf_counter()

print(
    "\nwith gprof2dot and graphviz (command dot):\n"
    "gprof2dot -f pstats profile.pstats | dot -Tpng -o profile.png"
)
