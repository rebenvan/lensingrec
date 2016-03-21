CC = gcc
OPTIMIZE = -g -O3

MLIBS = -lm

OBJS = main.o globle.o calculate.o readfile.o write_info.o
EXEC = lensingrec

all: $(EXEC)

lensingrec: $(OBJS)
	$(CC) $(OPTIMIZE) -o $@ $(OBJS) $(MLIBS)

$(OBJS): def.h Makefile
