.POSIX:
.SUFFIXES:
.SUFFIXES: .o .f90
include ../conf.mk
L = libncdff.a
all: $L
include obj.mk
include mod.mk
$L: $O; ar rv $@ $O && ranlib $@
.f90.o:; $(FC) -c $(FCFLAGS) $(FXFLAGS) $<
install: $L
	mkdir -p $(PREFIX)/lib $(PREFIX)/include && \
	cp $L $(PREFIX)/lib/ && \
	cp $M $(PREFIX)/include

clean:; rm -f -- $L $O $M
include dep.mk

# ./tool/list *.f90 *.inc | ./tool/dep > dep.mk
# ls *.f90 | awk 'BEGIN {print "O = \\"} {sub(/\.f90/, ".o"); print $0 "\\"} END {print ""}'  > obj.mk
# ls *.mod | awk 'BEGIN {print "M = \\"} {print $0 "\\"} END {print ""}'  > mod.mk
