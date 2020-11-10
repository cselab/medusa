.POSIX:
.SUFFIXES:
.SUFFIXES: .o .f90
include ../conf.mk
L = libppm.a
all: $L
# ls *.f90 | awk 'BEGIN {print "O = \\"} {sub(/\.f90/, ".o"); print $0 "\\"} END {print ""}'  > obj.mk
include obj.mk
0 = ${O:.o=.mod}
M = ${0:ppm_module_neighlist_MkNeighIdx.mod=ppm_module_neighlist_mkneighidx.mod}
$L: $O; ar rv $@ $O && ranlib $@
.f90.o:; $(FC) -c $(FCFLAGS) $(FXFLAGS) $(FFTW_FCFLAGS) $<
install: $L
	mkdir -p $(PREFIX)/lib $(PREFIX)/include && \
	cp $L $(PREFIX)/lib/ && \
	cp $M $(PREFIX)/include && \
	cp ppm_param.h $(PREFIX)/include

clean:; rm -f -- $L $O $M

# ./tool/list *.f90 *.h *.inc | ./tool/dep > dep.mk
include dep.mk
