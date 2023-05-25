all: build/plot.pdf

# hier Python-Skripte:
build/plot.pdf: build/A1 build/A2 | build
	python plot_A1.py
	python plot_A2.py

build/A1: A1.cpp | build
	g++ -O3 -o build/A1 A1.cpp
	./build/A1

build/A2: A2.cpp | build
	g++ -O3 -o build/A2 A2.cpp
	./build/A2

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean plots
