all: build/plot.pdf \
	build/A1_b_final.png \
	build/A1_c_final.png \
	build/A1_d_final.png \
	build/A1_e_final.png

# hier Python-Skripte:
build/plot.pdf: plot_A2.py build/A1 build/A2 | build
	python plot_A2.py

build/A1: A1.cpp | build
	g++ -O3 -o build/A1 A1.cpp
	./build/A1

build/A2: A2.cpp | build
	g++ -O3 -o build/A2 A2.cpp
	./build/A2

build/A1_%_final.png build/A1_%_anim.mp4: build/A1_%_E.txt build/A1_%_phi.txt
	python plot_A1.py $*

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean plots
