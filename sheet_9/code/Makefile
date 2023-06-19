all: build/A1_values.pdf

# hier Python-Skripte:
build/A1_values.pdf: plot.py build/A1 build/A2 | build
	python plot.py

build/A1: A1.cpp optimize.h | build
	g++ -O3 -o build/A1 A1.cpp
	./build/A1

build/A2: A2.cpp optimize.h | build
	g++ -O3 -o build/A2 A2.cpp
	./build/A2

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean plots
