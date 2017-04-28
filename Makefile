all: poisson_bin

.PHONY: all test clean format ready

model.hpp:
	curl "https://raw.githubusercontent.com/vlanore/tinycompo/master/model.hpp" > model.hpp

doctest.h:
	curl "https://raw.githubusercontent.com/onqtam/doctest/master/doctest/doctest.h" > doctest.h

poisson_bin: poissonGamma.cpp model.hpp doctest.h graphicalModel.hpp
	$(CXX) -std=gnu++11 $< -o $@ -DDOCTEST_CONFIG_DISABLE

test: poisson_bin
	@echo "============TEST=============" ; ./$< && echo "============================="

clean:
	rm -rf *_bin *.o model.hpp

format:
	@clang-format -i poissonGamma.cpp graphicalModel.hpp

ready: format test
	git status
