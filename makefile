SHELL = /bin/bash
MOVIE_FILE = trains.hvqm
MOVIE = samples/$(MOVIE_FILE)
ALL_MOVIES = \
	samples/*.hvqm \
	samples/pokemon/*.hvqm \
	samples/yakouchuu/*.hvqm

all: native

output_dir:
	mkdir -p /tmp/output
	rm -f output/*.ppm

build_emu:
	./toolchain/bin/mips-linux-gcc -Wall -Wextra -Og -g -I. hvqm2.c -L. -lhvqm2_pc -static -o hvqm2

emu: build_emu output_dir
	qemu-mips hvqm2 $(MOVIE)
	diff -qr output/ reference/$(MOVIE_FILE)/

debug: build_emu output_dir
	qemu-mips -g 1234 hvqm2 $(MOVIE) &
	toolchain/bin/mips-linux-gdb -ex 'target remote localhost:1234' -ex c hvqm2

build_native:
	#clang -DNATIVE=1 -Wall -Wextra -Og -g -I. hvqm2.c -o hvqm2 -fsanitize=address
	gcc -DNATIVE=1 -Wall -Wextra -Og -g -I. hvqm2.c -o hvqm2 -fsanitize=address
	#gdb -ex r --args ./hvqm2 $(MOVIE)

run: build_native output_dir
	./hvqm2 $(MOVIE)
	diff -qr output/ reference/$(MOVIE_FILE)/

test: build_native output_dir
	for i in $(ALL_MOVIES); do echo $$i; rm -f output/*.ppm; ./hvqm2 $$i && diff -qr output/ reference/$${i:8}/; done
