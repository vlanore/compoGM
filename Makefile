TINYCOMPO_FILES = tinycompo.hpp arrays.hpp

all: poisson_bin sample_bin

.PHONY: all test clean format ready

$(TINYCOMPO_FILES):
	curl "https://raw.githubusercontent.com/vlanore/tinycompo/master/tinycompo.hpp" > tinycompo.hpp
	curl "https://raw.githubusercontent.com/vlanore/tinycompo/master/arrays.hpp" > arrays.hpp

poisson_bin: poissonGamma.cpp $(TINYCOMPO_FILES) graphicalModel.hpp
	$(CXX) -std=gnu++11 $< -o $@

sample_bin: samplePlate.cpp $(TINYCOMPO_FILES) graphicalModel.hpp
	$(CXX) -std=gnu++11 $< -o $@

test: poisson_bin
	@echo "============TEST=============" ; ./$< && echo "============================="

clean:
	rm -rf *_bin *.o $(TINYCOMPO_FILES)

format:
	@clang-format -i poissonGamma.cpp graphicalModel.hpp

ready: format test
	git status
