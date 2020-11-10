CC=gcc
CFLAGS=-gdwarf-2 -g
OBJS=file-io.o load_i2.o

all: main

main: $(OBJS)
	@echo Building main...
	$(CC) $(CFLAGS) $(OBJS) -lz test_main.c -o main
	@echo Done.

file-io.o: file-io.h file-io.c
load_i2.o: load_i2.h load_i2.c

pileup.o : pileup.h pileup.c
	echo "Making pileup.o..."
	$(CC) $(CFLAGS) -c pileup.c -lz -o pileup.o

file-io.o : file-io.h file-io.c
	echo "Making file-io.o..."
	$(CC) $(CFLAGS) -c file-io.c -lz -o file-io.o
