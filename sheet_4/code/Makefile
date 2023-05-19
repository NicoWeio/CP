all: build/plot.pdf \
	build/A1_b_final.pdf \
	build/A1_c_final.pdf \
	build/A1_d_final.pdf \
	build/A1_e_final.pdf \
	build/A1_c_analytical.pdf \

# hier Python-Skripte:
build/plot.pdf: plot_A2.py build/A1 build/A2 | build
	python plot_A2.py

build/A1_c_analytical.pdf: A1c.py | build
	python A1c.py

build/A1: A1.cpp | build
	g++ -O3 -o build/A1 A1.cpp
	./build/A1

build/A2: A2.cpp | build
	g++ -O3 -o build/A2 A2.cpp
	./build/A2

build/A1_%_final.pdf build/A1_%_anim.mp4: plot_A1.py build/A1_%_E.txt build/A1_%_phi.txt | build
	python plot_A1.py $*

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean plots
