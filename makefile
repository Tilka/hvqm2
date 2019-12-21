MOVIE_FILE = sample.hvqm
MOVIE = samples/$(MOVIE_FILE)

all: native

build_emu:
	./toolchain/bin/mips-linux-gcc -Wall -Wextra -Og -g -I. hvqm2.c -L. -lhvqm2_pc -static -o hvqm2
	mkdir -p /tmp/output
	rm -f output/*.ppm

emu: build_emu
	qemu-mips hvqm2 $(MOVIE)
	diff -qr output/ reference/$(MOVIE_FILE)/

debug: build_emu
	qemu-mips -g 1234 hvqm2 $(MOVIE) &
	toolchain/bin/mips-linux-gdb -ex 'target remote localhost:1234' -ex c hvqm2

native:
	#clang -DNATIVE=1 -Wall -Wextra -Og -g -I. hvqm2.c -o hvqm2 -fsanitize=address
	gcc -DNATIVE=1 -Wall -Wextra -Og -g -I. hvqm2.c -o hvqm2 -fsanitize=address
	#gdb -ex r --args ./hvqm2 $(MOVIE)
	mkdir -p /tmp/output
	rm -f output/*.ppm

run: native
	./hvqm2 $(MOVIE)
	diff -qr output/ reference/$(MOVIE_FILE)/
