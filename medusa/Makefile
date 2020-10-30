.POSIX:
.SUFFIXES: .o .f
L = ppm_tv
all: $L
# ls ppm_module_*.f | awk 'BEGIN {print "O = \\"} {sub(/\.f/, ".o"); print $0 "\\"} END {print ""}'  > obj.mk
include obj.mk

PREFIX = $(HOME)/.local
FCFLAGS = -O2 -g
PPM_FCFLAGS = -I$(PREFIX)/include
FC = mpifort
LINK = $(FC)
M = ${O:.o=.mod}
M_FCFLAGS = -cpp -ffree-line-length-none -fallow-argument-mismatch -ffree-form
M_LDFLAGS = -L$(PREFIX)/lib -Wl,-R$(PREFIX)/lib -lppm -lnetcdff -lfftw3f
$L: $O; $(LINK) $O $(LDFLAGS) $(M_LDFLAGS) -o $@
.f.o:; $(FC) -c $(FCFLAGS) $(M_FCFLAGS) $(PPM_FCFLAGS) -c $<
clean:; rm -f -- $L $O $M

# ./tool/list *.f *.h *.inc | ./tool/dep > dep.mk
include dep.mk
