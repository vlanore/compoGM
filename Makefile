.PHONY: all clean ready test

all: test_bin

tinycompo.hpp:
	curl https://raw.githubusercontent.com/vlanore/tinycompo/master/tinycompo.hpp > $@

%_bin: src/%.cpp tinycompo.hpp
	$(CXX) --std=c++11 $< -o $@

clean:
	rm *_bin

test: test_bin
	./$<

format:
	clang-format -i src/*.hpp src/*.cpp

ready: all
	make test --no-print-directory
	@git status
