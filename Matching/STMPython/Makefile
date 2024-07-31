PATH_INPUT_DATA="../../Documentation/TestData/Processed_DATA/MyExperiment/Parallel/Matching/Rays/rays_11-20.dat"

all:
	transonic stm_util.py -af "-march=native -Ofast -DUSE_XSIMD"
	transonic util_groupby.py -af "-march=native -Ofast -DUSE_XSIMD"

clean:
	rm -rf __pythran__

bench:
	./stm_hdf5.py $(PATH_INPUT_DATA) 1 2 2 0.2 400 400 250 2

bench_old:
	STM_PYTHON_USE_UNOPTIMIZED=1 ./stm_hdf5.py $(PATH_INPUT_DATA) 1 2 2 0.2 400 400 250 2

bench_cpp:
	../STMCpp/STM -i $(PATH_INPUT_DATA) -f 1 -c 2 -d 0.2 -m 2 -x 400 -y 400 -z 250 -b -140 140 -150 150 5 170 --hdf5

bench_min-distance-matches-1ray:
	./stm_hdf5.py $(PATH_INPUT_DATA) 1 2 2 0.2 400 400 250 2 -md1r 2

profile:
	python profile_simple.py

black:
	black -l 82 .

test:
	pytest test_stm_util.py
