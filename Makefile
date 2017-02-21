all: test_bin

test_bin: test.cpp model.hpp
	$(CXX) -std=gnu++11 $< -o test_bin

test: test_bin
	./test_bin

clean:
	rm -rf *_bin *.o

format:
	clang-format -i test.cpp model.hpp
