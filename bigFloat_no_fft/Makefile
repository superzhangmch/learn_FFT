CXX = clang++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

all: calc math_test pi1024

calc: calc.cpp bfmath.h bigfloat.h
	$(CXX) $(CXXFLAGS) -o $@ calc.cpp

math_test: math_test.cpp bfmath.h bigfloat.h
	$(CXX) $(CXXFLAGS) -o $@ math_test.cpp

pi1024: pi1024.cpp bigfloat.h
	$(CXX) $(CXXFLAGS) -o $@ pi1024.cpp

clean:
	rm -f calc math_test pi1024

.PHONY: all clean
