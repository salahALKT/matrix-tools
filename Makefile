CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

DEPS = Matrix.h Matrix.tpp overloads.h

.PHONY: all test clean

all: demo

demo: main.cpp overloads.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) main.cpp overloads.cpp -o $@

tests/test_matrix: tests/test_matrix.cpp overloads.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) tests/test_matrix.cpp overloads.cpp -o $@

test: tests/test_matrix
	./tests/test_matrix

clean:
	rm -f demo tests/test_matrix *.o
