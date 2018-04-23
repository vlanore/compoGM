CPPFLAGS= -Wall -Wextra -O3 --std=c++11

.PHONY: all clean ready test format

all: test_bin M1_bin

tinycompo.hpp:
	@echo "-- Downloading tinycompo.hpp from github..."
	curl https://raw.githubusercontent.com/vlanore/tinycompo/experimental/tinycompo.hpp > $@
	@echo "-- Done."

%_bin: src/%.cpp tinycompo.hpp src/*.hpp
	$(CXX) -I. $(CPPFLAGS) $< -o $@

clean:
	rm -f *_bin *.hpp

test: test_bin
	./$<

m1: M1_bin
	./$<

format:
	clang-format -i src/*.hpp src/*.cpp

ready:
	@echo "-- Formatting with clang-format..."
	@make format --no-print-directory
	@echo "\n-- Compiling if necessary..."
	@make --no-print-directory
	@echo "\n-- Launching test..."
	@make test --no-print-directory
	@echo "\n-- All done, git status is:"
	@git status
