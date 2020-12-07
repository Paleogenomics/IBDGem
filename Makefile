CC=gcc
CFLAGS=-gdwarf-2 -g
OBJS=file-io.o load_i2.o pileup.o nchoosek.o

all: main

main: $(OBJS)
	@echo Building main...
	$(CC) $(CFLAGS) $(OBJS) -lz compare-vcf2bam.c -o compare-vcf2bam
	@echo Done.

file-io.o: file-io.h file-io.c
load_i2.o: load_i2.h load_i2.c
pileup.o: pileup.h pileup.c
nchoosek.o: nchoosek.h nchoosek.c
