all: test_bin poisson_bin
mpi: test_mpi_bin

.PHONY: all mpi test testmpi clean

test_bin: test.cpp model.hpp
	$(CXX) -std=gnu++11 $< -o $@

poisson_bin: poissonGamma.cpp model.hpp
	$(CXX) -std=gnu++11 $< -o $@

test_mpi_bin: test.cpp model.hpp
	mpic++ -std=gnu++11 $< -o test_mpi_bin -DUSE_MPI

test: poisson_bin
	./$<

testmpi: test_mpi_bin
	mpirun -np 2 ./$<

dot: tmp.dot
	@dot -Tps $< -o tmp.ps
	@evince tmp.ps &

clean:
	rm -rf *_bin *.o

format:
	clang-format -i test.cpp model.hpp
