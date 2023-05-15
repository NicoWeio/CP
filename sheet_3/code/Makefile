all: build/A1

# hier Python-Skripte:
#build/plot.pdf: plot.py build/A1 | build
#	python plot.py

build/A1: A1.cpp | build
	g++ -O3 -o build/A1 A1.cpp
	./build/A1

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean plots
