CPPFLAGS= -Wall -Wextra -Wfatal-errors -O3 --std=c++11 -pthread -march=native

all: test_bin M0_bin M0_mpi_bin M1_bin M2_bin M3_bin

tinycompo.hpp:
	@echo "Downloading tinycompo.hpp from github..."
	curl https://raw.githubusercontent.com/vlanore/tinycompo/experimental/tinycompo.hpp > $@
	@echo "Done."

csv-parser.hpp:
	@echo "Downloading csv parser from github..."
	curl https://raw.githubusercontent.com/AriaFallah/csv-parser/master/parser.hpp > $@
	@echo "Done."

%_bin: src/%.cpp tinycompo.hpp csv-parser.hpp src/*.hpp
	mpic++ -I. $(CPPFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f *_bin *.hpp

.PHONY: test
test: test_bin
	./$<

.PHONY: m0
m0: M0_bin
	./$<

.PHONY: m0_mpi
m0_mpi: M0_mpi_bin
	mpirun -np 3 ./$<

.PHONY: m1
m1: M1_bin
	./$< ~/data/rnaseq_mini

.PHONY: m2
m2: M2_bin
	./$< ~/data/rnaseq_mini

.PHONY: m3
m3: M3_bin
	./$< ~/data/rnaseq_mini

.PHONY: tmp
tmp: tmp_bin
	./$<

.PHONY: format
format:
	clang-format -i src/*.hpp src/*.cpp

.PHONY: ready
ready:
	@echo "\033[1m\033[95mFormatting with clang-format...\033[0m"
	@make format --no-print-directory
	@echo "\033[1m\033[95m\nCompiling if necessary...\033[0m"
	@make -j6 --no-print-directory
	@echo "\033[1m\033[95m\nLaunching test...\033[0m"
	@make test --no-print-directory
	@echo "\033[1m\033[95m\nAll done, git status is:\033[0m"
	@git status
