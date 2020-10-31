.POSIX:
.SUFFIXES: .o .c
include ../conf.mk
L = libmetis.a
all: $L
include obj.mk

PREFIX = $(HOME)/.local
LINK = $(CC)
$L: $O; ar rv $@ $O && ranlib $@
.c.o:; $(CC) -c $(CFLAGS) -I. $<
install: $L
	mkdir -p $(PREFIX)/lib && \
	cp $L $(PREFIX)/lib

clean:; rm -f -- $L $O

include dep.mk
