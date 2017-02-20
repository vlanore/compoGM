all: test_bin

test_bin: test.cpp
	$(CXX) $< -o test_bin

test: test_bin
	./test_bin

clean:
	rm -rf *_bin *.o
