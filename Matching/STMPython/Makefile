PATH_INPUT_DATA="../../Documentation/TestData/Processed_DATA/MyExperiment/Parallel/Matching/Rays/rays_11-20.dat"

install: sync build

build:
	pdm run transonic stm_util.py -af "-march=native -Ofast -DUSE_XSIMD"
	pdm run transonic util_groupby.py -af "-march=native -Ofast -DUSE_XSIMD"

clean:
	rm -rf __pythran__

bench:
	pdm run ./stm_hdf5.py $(PATH_INPUT_DATA) 1 2 2 0.2 400 400 250 2

bench_unoptimized:
	STM_PYTHON_USE_UNOPTIMIZED=1 pdm run ./stm_hdf5.py $(PATH_INPUT_DATA) 1 2 2 0.2 400 400 250 2

bench_cpp:
	../STMCpp/STM -i $(PATH_INPUT_DATA) -f 1 -c 2 -d 0.2 -m 2 -x 400 -y 400 -z 250 -b -140 140 -150 150 5 170 --hdf5

bench_min-distance-matches-1ray:
	pdm run ./stm_hdf5.py $(PATH_INPUT_DATA) 1 2 2 0.2 400 400 250 2 -md1r 2

profile:
	pdm run python profile_simple.py

black:
	black -l 82 .

test:
	pdm run pytest test_stm_util.py

lock:
	pdm lock --exclude-newer 2022-08-01

sync:
	pdm sync
