PATH_INPUT_DATA="../../Documentation/TestData/Processed_DATA/MyExperiment/Parallel/Matching/Rays/rays_1-10.dat"

all:
	transonic stm_util.py

clean:
	rm -rf __pythran__

bench:
	python stm.py $(PATH_INPUT_DATA) 1 2 2 0.2 400 400 250 2

bench_old:
	STM_PYTHON_USE_UNOPTIMIZED=1 python stm.py $(PATH_INPUT_DATA) 1 2 2 0.2 400 400 250 2

bench_cpp:
	../STMCpp/STM -i $(PATH_INPUT_DATA) -f 1 -c 2 -d 0.2 -m 2 -x 400 -y 400 -z 250 -b -140 140 -150 150 5 170 --hdf5

profile:
	python profile_simple.py