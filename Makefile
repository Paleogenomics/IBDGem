CC=gcc
CFLAGS=-ggdb3 -Wall
OBJS=$(addprefix $(BUILD)/, file-io.o pileup.o ibd-parse.o ibd-math.o)
LDFLAGS=-lz -lm

SRC   := src
BUILD := build
vpath %.c $(SRC)

.PHONY: all
all: ibdgem hiddengem
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

.PHONY: clean
clean:
	@echo Removing all object files and executables...
	rm -f $(OBJS) ibdgem hiddengem
	@echo Done.