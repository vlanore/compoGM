CPPFLAGS= -Wall -Wextra -O3 --std=c++11 -pthread

.PHONY: all clean ready test format tmp m1 m2

all: test_bin M1_bin M2_bin

tinycompo.hpp:
	@echo "-- Downloading tinycompo.hpp from github..."
	curl https://raw.githubusercontent.com/vlanore/tinycompo/experimental/tinycompo.hpp > $@
	@echo "-- Done."

csv-parser.hpp:
	@echo "-- Downloading csv parser from github..."
	curl https://raw.githubusercontent.com/AriaFallah/csv-parser/master/parser.hpp > $@
	@echo "-- Done."

%_bin: src/%.cpp tinycompo.hpp csv-parser.hpp src/*.hpp
	$(CXX) -I. $(CPPFLAGS) $< -o $@

clean:
	rm -f *_bin *.hpp

test: test_bin
	./$<

m1: M1_bin
	./$< ~/data/rnaseq_mini

m2: M2_bin
	./$< ~/data/rnaseq_mini

m3: M3_bin
	./$< ~/data/rnaseq_mini

tmp: tmp_bin
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
