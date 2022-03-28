CC=gcc
CFLAGS=-g -Wall -pthread
OBJS=file-io.o load_i2.o pileup.o nchoosek.o
LDFLAGS=-lz -lm

.PHONY: all
all: ibdgem hiddengem aggregate

ibdgem: $(OBJS)
	@echo Building ibdgem...
	$(CC) $(CFLAGS) $(OBJS) $@.c $(LDFLAGS) -o $@
	@echo Done.

hiddengem: file-io.o
	@echo Building hiddengem...
	$(CC) $(CFLAGS) $? $@.c -lz -o $@
	@echo Done.

aggregate: file-io.o
	@echo Building aggregate...
	$(CC) $(CFLAGS) $? $@.c -lz -o $@
	@echo Done.

file-io.o: file-io.h file-io.c
load_i2.o: load_i2.h load_i2.c
pileup.o: pileup.h pileup.c
nchoosek.o: nchoosek.h nchoosek.c

.PHONY: clean
clean:
	@echo Removing all object files and executables...
	rm -f $(OBJS) ibdgem hiddengem aggregate
	@echo Done.