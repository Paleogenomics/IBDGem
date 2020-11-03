CC=gcc
CFLAGS=-gdwarf-2 -g
OBJS=file-io.o load_i2.o

all: main

main: $(OBJS)
	@echo Building main...
	$(CC) $(CFLAGS) test_main.c -o main -lz $(OBJS)
	@echo Done.

file-io.o: file-io.h file-io.c
load_i2.o: load_i2.h load_i2.c