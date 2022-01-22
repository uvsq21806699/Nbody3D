# -mavx2 don't work on my computer'

all: nbody.g nbody.i

nbody_fast.g: nbody.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp 

nbody_0.g: nbody.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp 

nbodyv.g: nbody.c
	gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp -DVERBOSE

nbody.i: nbody.c
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

nbodyv.i: nbody.c
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp -DVERBOSE

local: nbody.c
	gcc -march=native -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o nbody.g -lm -fopenmp
	gcc -march=native -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o nbodyv.g -lm -fopenmp -DVERBOSE

clean:
	rm -Rf *~ nbody.g nbody_fast.g nbody_0.g nbody.i nbodyv.g nbodyv.i *.optrpt Final_pos2.txt