CC=gcc
CFLAGS=-gdwarf-2 -g

file-io.o : file-io.h file-io.c
	echo "Making file-io..."
	$(CC) $(CFLAGS) -c file-io.c -lz -o file-io.o
