all: build/plot.pdf build/A1_r.pdf

build/A1_r.pdf: plot.py build/A1 build/A2 | build
	./build/A1 && python plot.py

# here are the Python scripts:
build/plot.pdf: A2_plot.py build/A2 | build
	./build/A2 && python A2_plot.py

# here are the C++ programs:
build/A1: A1.cpp | build
	g++ -o build/A1 A1.cpp

build/A2: A2.cpp | build
	g++ -o build/A2 A2.cpp


# create the build directory if it doesn't exist
build:
	mkdir -p build

# clean up the build directory
clean:
	rm -rf build

.PHONY: all clean
