all: test_bin
mpi: test_mpi_bin

.PHONY: all mpi test testmpi clean

test_bin: test.cpp model.hpp
	mpic++ -std=gnu++11 $< -o test_bin

test_mpi_bin: test.cpp model.hpp
	mpic++ -std=gnu++11 $< -o test_mpi_bin -DUSE_MPI

test: test_bin
	./$<

testmpi: test_mpi_bin
	mpirun -np 2 ./$<

clean:
	rm -rf *_bin *.o

format:
	clang-format -i test.cpp model.hpp
