all: build/plot.pdf

# hier Python-Skripte:
build/plot.pdf: plot.py build/A1 build/A2 build/A3 | build
	python plot.py

build/A1: A1.cpp | build
	g++ -O3 -o build/A1 A1.cpp
	./build/A1

build/A2: A2.cpp | build
	g++ -O3 -o build/A2 A2.cpp
	./build/A2

build/A3: A3.cpp | build
	g++ -O3 -o build/A3 A3.cpp
	./build/A3

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean plots
