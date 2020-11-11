.POSIX:
.SUFFIXES:
.SUFFIXES: .f90
include ../conf.mk

M = main
all: $M
.f90:
	$(FC) $< $(FCFLAGS) $(FXFLAGS) -I$(PREFIX)/include -L$(PREFIX)/lib -lfftw3 -lfftw3f -o $@
clean:; rm -f $M
