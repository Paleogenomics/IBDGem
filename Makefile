CC=gcc
CFLAGS=-ggdb3 -Wall -pthread
OBJS=$(addprefix $(BUILD)/, file-io.o load-i2.o pileup.o ibd-math.o)
LDFLAGS=-lz -lm

SRC   := src
BUILD := build
vpath %.c $(SRC)

.PHONY: all
all: ibdgem hiddengem aggregate

$(BUILD)/%.o: %.c
	mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

ibdgem: $(OBJS)
	@echo Building ibdgem...
	$(CC) $(CFLAGS) $(OBJS) $(SRC)/$@.c $(LDFLAGS) -o $@
	@echo Done.

hiddengem: $(BUILD)/file-io.o
	@echo Building hiddengem...
	$(CC) $(CFLAGS) $? $(SRC)/$@.c -lz -o $@
	@echo Done.

aggregate: $(BUILD)/file-io.o
	@echo Building aggregate...
	$(CC) $(CFLAGS) $? $(SRC)/$@.c -lz -o $@
	@echo Done.

.PHONY: clean
clean:
	@echo Removing all object files and executables...
	rm -f $(OBJS) ibdgem hiddengem aggregate
	@echo Done.
