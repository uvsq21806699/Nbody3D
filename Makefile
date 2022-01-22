all:
	build

name = t

build/test: v1/test.txt
	name = $<
	sed -i 's|/|_|g' $(name)
#	echo $(name)

	echo $@ $< $? $* $(dir $@)$(name)
	cp -f $< $(dir $@)$(sed 's|/|_|g' $<)

build/nbody.g: 0/nbody.c
	echo "$@ $< $? $*"
#gcc -march=native -mavx2 -Ofast -fopt-info-all=nbody.gcc.optrpt $< -o $@ -lm -fopenmp

build/nbody.i: 0/nbody.c
	icc -xhost -Ofast -qopt-report $< -o $(basename src/foo.c src-1.0/bar hacks) -qmkl -qopenmp

build/nbody.i: v1/nbody.c
	icc -xhost -Ofast -qopt-report $< -o $@ -qmkl -qopenmp

dir:
	mkdir -p build


performance:
	sudo cpupower -c all frequency-set -g performance

powersave:
	sudo cpupower -c all frequency-set -g powersave

