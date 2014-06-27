CC = gcc
CFLAGS = -pthread -std=gnu99 -Og -g -Wall -Wextra -pedantic -m64

TESTS = fmitest
PROGS = search_reads build_index

all: $(TESTS) $(PROGS)

build_index: seqindex.o csacak.o build_index.o fileio.o

search_reads: seqindex.o csacak.o search_reads.o fileio.o
	gcc -o $@ $^ $(CFLAGS)

fmitest: fmitest.o seqindex.o csacak.o
	gcc -o $@ $^ $(CFLAGS)

# No, gcc, I will not listen to your whinging
csacak.o: csacak.c
	gcc -std=gnu99 -O3 -m64 -c $^

clean:
	rm -f *.o $(TESTS) $(PROGS) *~

.PHONY: clean all
