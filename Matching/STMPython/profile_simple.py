"""
Example of output:

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1   41.375   41.375   44.827   44.827 stm_util.py:574(<listcomp>)
 26944128   36.134    0.000   42.883    0.000 stm_util.py:21(<genexpr>)
   549354   22.618    0.000   27.463    0.000 {built-in method builtins.sorted}
        1   12.385   12.385  231.874  231.874 stm_util.py:375(space_traversal_matching)
    10429   10.003    0.001   17.754    0.002 stm_util.py:59(uniquify)
    10428    8.134    0.001    8.134    0.001 stm_util.py:566(<listcomp>)
  2568545    7.490    0.000    7.984    0.000 stm_util.py:583(<listcomp>)
 70728336    6.748    0.000    6.748    0.000 stm_util.py:21(<lambda>)
    61386    5.820    0.000    5.820    0.000 stm_util.py:214(<listcomp>)
  1813940    5.207    0.000    9.784    0.000 stm_util.py:595(<lambda>)
 30329636    4.830    0.000    4.830    0.000 stm_util.py:570(<lambda>)
  3368016    4.829    0.000   47.712    0.000 stm_util.py:19(expand_neighbours)

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
