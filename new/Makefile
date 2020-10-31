.POSIX:
.SUFFIXES:
.SUFFIXES: .o .f90
include ../conf.mk
L = libncdf_cse.a
all: $L
include obj.mk
include mod.mk
#M = ${O:.o=.mod}
$L: $O; ar rv $@ $O && ranlib $@
.f90.o:; $(FC) -c $(FCFLAGS) $(FXFLAGS) $<
install: $L
	mkdir -p $(PREFIX)/lib $(PREFIX)/include && \
	cp $L $(PREFIX)/lib/ && \
	cp $M $(PREFIX)/include

clean:; rm -f -- $L $O $M

# ./tool/list *.f90 *.h *.inc | ./tool/dep > dep.mk
#include dep.mk
# ls *.f90 | awk 'BEGIN {print "O = \\"} {sub(/\.f90/, ".o"); print $0 "\\"} END {print ""}'  > obj.mk
# ls *.mod | awk 'BEGIN {print "M = \\"} {print $0 "\\"} END {print ""}'  > mod.mk
