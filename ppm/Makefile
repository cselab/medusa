.POSIX:
.SUFFIXES: .o .f
L = libppm.a
all: $L
# ls ppm_module_*.f | awk 'BEGIN {print "O = \\"} {sub(/\.f/, ".o"); print $0 "\\"} END {print ""}'  > obj.mk
include obj.mk

PREFIX = $(HOME)/.local
FCFLAGS = -O2 -g
FC = mpifort
0 = ${O:.o=.mod}
M = ${0:ppm_module_neighlist_MkNeighIdx.mod=ppm_module_neighlist_mkneighidx.mod}

M_FCFLAGS = -cpp -ffree-line-length-none -fallow-argument-mismatch -ffree-form
$L: $O; ar rv $@ $O && ranlib $@
.f.o:; $(FC) -c $(FCFLAGS) $(M_FCFLAGS) $<
install: $L
	mkdir -p $(PREFIX)/lib $(PREFIX)/include && \
	cp $L $(PREFIX)/lib/ && \
	cp $M $(PREFIX)/include && \
	cp ppm_param.h $(PREFIX)/include

clean:; rm -f -- $L $O $M

# ./tool/list *.f *.h *.inc | ./tool/dep > dep.mk
include dep.mk
