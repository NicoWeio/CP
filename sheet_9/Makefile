all: build/main.pdf

# hier Python-Skripte:
plot: FORCE
	cd code && make

# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/main.pdf: plot

build/main.pdf: FORCE | build
	  TEXINPUTS=build: \
#	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	main.tex
	evince build/main.pdf &

build:
	mkdir -p build

clean:
	rm -rf build
	cd code && make clean

FORCE:

.PHONY: all clean
