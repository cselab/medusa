.POSIX:
.SUFFIXES:
.SUFFIXES: .o .f90
include ../conf.mk
L = ppm_tv
all: $L
# ls *.f90 | awk 'BEGIN {print "O = \\"} {sub(/\.f90/, ".o"); print $0 "\\"} END {print ""}'  > obj.mk
include obj.mk

PPM_FCFLAGS = -I$(PREFIX)/include
PPM_LDFLAGS = -L$(PREFIX)/lib -lppm
LINK = $(FC)
M = module_wvic.mod ppm_module_fft.mod wvic_module_io.mod
M_LDFLAGS = -lnetcdff -lnetcdf
$L: $O; $(LINK) $O $(LDFLAGS) $(PPM_LDFLAGS) $(FFTW_LDFLAGS) $(M_LDFLAGS) -o $@
.f90.o:; $(FC) $< -c $(FCFLAGS) $(FXFLAGS) $(PPM_FCFLAGS) $(FFTW_FCFLAGS)
clean:; rm -f -- $L $O $M

# ./tool/list *.f90 *.inc | ./tool/dep > dep.mk
include dep.mk
