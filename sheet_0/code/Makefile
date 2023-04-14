all: build/a3_a.pdf

# define a variable for g++ flags
CXXFLAGS = -Wall -Wextra -Werror -pedantic -std=c++11 -O3

# hier Python-Skripte:
build/a3_a.pdf: plot.py build/A1_hello_world build/A2_rounding_error build/A3 | build
	python plot.py

build/A1_hello_world: A1_hello_world.cpp | build
	g++ $(CXXFLAGS) -o build/A1_hello_world A1_hello_world.cpp
	./build/A1_hello_world

build/A2_rounding_error: A2_rounding_error.cpp | build
	g++ $(CXXFLAGS) -o build/A2_rounding_error A2_rounding_error.cpp
	./build/A2_rounding_error

build/A3: A3.cpp | build
	g++ $(CXXFLAGS) -o build/A3 A3.cpp
	./build/A3

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean plots
