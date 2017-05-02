all: poisson_bin sample_bin

.PHONY: all test clean format ready

model.hpp:
	curl "https://raw.githubusercontent.com/vlanore/tinycompo/master/model.hpp" > model.hpp

poisson_bin: poissonGamma.cpp model.hpp doctest.h graphicalModel.hpp
	$(CXX) -std=gnu++11 $< -o $@ -DDOCTEST_CONFIG_DISABLE

sample_bin: samplePlate.cpp model.hpp doctest.h graphicalModel.hpp
	$(CXX) -std=gnu++11 $< -o $@ -DDOCTEST_CONFIG_DISABLE

test: poisson_bin
	@echo "============TEST=============" ; ./$< && echo "============================="

clean:
	rm -rf *_bin *.o model.hpp

format:
	@clang-format -i poissonGamma.cpp graphicalModel.hpp

ready: format test
	git status
