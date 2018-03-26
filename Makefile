CPPFLAGS= -Wall -Wextra -O3 --std=c++11

.PHONY: all clean ready test

all: test_bin

tinycompo.hpp:
	curl https://raw.githubusercontent.com/vlanore/tinycompo/experimental/tinycompo.hpp > $@

tinycompo_mpi.hpp:
	curl https://raw.githubusercontent.com/vlanore/tinycompo/experimental/tinycompo.hpp > $@

%_bin: src/%.cpp tinycompo.hpp tinycompo_mpi.hpp src/interfaces.hpp
	$(CXX) -I. $(CPPFLAGS) $< -o $@

clean:
	rm -f *_bin tinycompo.hpp

test: test_bin
	./$<

format:
	clang-format -i src/*.hpp src/*.cpp

ready: all
	make test --no-print-directory
	@git status
