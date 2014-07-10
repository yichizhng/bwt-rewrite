CC = gcc
CFLAGS = -pthread -std=gnu99 -O3 -m64

TESTS = fmitest filetest
PROGS = search_reads build_index single_align

all: $(TESTS) $(PROGS)

single_align: single_align.c csacak.o fileio.o seqindex.o smw.o stack.o
	gcc -o $@ $^ $(CFLAGS)

build_index: seqindex.o csacak.o build_index.o fileio.o
	gcc -o $@ $^ $(CFLAGS)

search_reads: seqindex.o csacak.o search_reads.o fileio.o stack.o
	gcc -o $@ $^ $(CFLAGS)

fmitest: fmitest.o seqindex.o csacak.o
	gcc -o $@ $^ $(CFLAGS)

filetest: filetest.o seqindex.o csacak.o fileio.o
	gcc -o $@ $^ $(CFLAGS)

# No, gcc, I will not listen to your whinging
csacak.o: csacak.c
	gcc -std=gnu99 -O3 -m64 -c $^

clean:
	rm -f *.o $(TESTS) $(PROGS) *~

.PHONY: clean all
