CC = gcc
OPTIMIZE = -g -O3

MLIBS = -lm
FITS = -lcfitsio
GSL = -lgsl -lgslcblas

FITS_LIB = -L/home/FAN/src/cfitsio-3.380/lib
FITS_INC = -I/home/FAN/src/cfitsio-3.380/include
GSL_LIB = -L/home/FAN/src/gsl-1.16/lib
GSL_INC = -I/home/FAN/src/gsl-1.16/include

CFLAGS = $(OPTIMIZE) $(FITS_INC) $(GSL_INC)
CLIBS = $(FITS_LIB) $(GSL_LIB) $(FITS) $(GSL)

OBJS1 = makegal.o healpix_utils.o
OBJS2 = makedeltat.o healpix_utils.o
EXEC = makegal makedeltat

all: $(EXEC)

makegal: $(OBJS1)
	$(CC) $(CFLAGS) -o $@ $(OBJS1) $(CLIBS) $(MLIBS)
makedeltat: $(OBJS2)
	$(CC) $(OPTIMIZE) $(GSL_INC) -o $@ $(OBJS2) $(GSL_LIB) $(GSL) $(MLIBS)

$(OBJS1): healpix_utils.h Makefile
$(OBJS2): healpix_utils.h Makefile
