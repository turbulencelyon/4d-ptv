MP_H5 = /opt/local/bin/h5c++

ifneq ("$(wildcard $(MP_H5))", "")
	CXX = $(MP_H5)
else
	CXX = /usr/bin/h5c++
endif

HEADERS = $(wildcard *.h)
SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))
RAYS = $(wildcard ../../Documentation/TestData/Processed_DATA/MyExperiment/Parallel/Matching/Rays/*.dat)
HDF5 = $(patsubst %.dat, %_out_cpp.h5, $(RAYS))
BIN = $(patsubst %.dat, %_out_cpp.bin, $(RAYS))
LOG = $(patsubst %.dat, %.log, $(RAYS))

CXXFLAGS = -std=gnu++11 -O3

all: STM STM_info

test: $(LOG)

%.o: %.cpp $(HEADERS)
	$(CXX) -c $(CXXFLAGS) $< -o $@

%_out_cpp.h5: %.dat STM
	./STM --config test.ini -i $<

%.log: %_out_cpp.h5 STM_info
	./STM_info out/`basename $<` > out/`basename $@`
	python STM_reader.py out/`basename $<` >> out/`basename $@`
	tail -n 1 $@
	
STM_info: STM_info.o STM_hdf5.o
	$(CXX) -o STM_info $(CXXFLAGS) STM_info.o STM_hdf5.o

STM: STM.o STM_hdf5.o STM_main.o STM_CMatrix.o STM_helpers.o
	$(CXX) -o STM $(CXXFLAGS) STM.o STM_hdf5.o STM_main.o STM_CMatrix.o STM_helpers.o

clean:
	rm -f STM STM_info $(OBJECTS)
	rm -f $(HDF5) $(BIN)

testlocal: STM STM_info
	./STM --config rays.ini -i ../Benchmark/rays_1-10.dat;
	python STM_reader.py rays_1-10;
	./STM --config rays.ini -i ../Benchmark/rays_21-30.dat;
	python STM_reader.py rays_21-30;
	./STM --config rays.ini -i ../Benchmark/rays_31-40.dat
	./STM_info rays_31-40_out_cpp.h5
	python STM_reader.py rays_31-40

