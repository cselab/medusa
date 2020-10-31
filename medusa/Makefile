.POSIX:
.SUFFIXES:
.SUFFIXES: .o .f90
include ../conf.mk
L = ppm_tv
all: $L
# ls *.f90 | awk 'BEGIN {print "O = \\"} {sub(/\.f90/, ".o"); print $0 "\\"} END {print ""}'  > obj.mk
include obj.mk

PPM_FCFLAGS = -I$(PREFIX)/include
PPM_LDFLAGS = -L$(PREFIX)/lib -Wl,-R$(PREFIX)/lib -lppm
LINK = $(FC)
M = ${O:.o=.mod}
M_LDFLAGS = -lnetcdff -lfftw3f
$L: $O; $(LINK) $O $(PPM_LDFLAGS) $(LDFLAGS) $(M_LDFLAGS) -o $@
.f90.o:; $(FC) -c $(FCFLAGS) $(FXFLAGS) $(PPM_FCFLAGS) -c $<
clean:; rm -f -- $L $O $M

# ./tool/list *.f90 *.inc | ./tool/dep > dep.mk
include dep.mk
