.POSIX:
.SUFFIXES: .o .c
L = libmetis.a
all: $L
include obj.mk

PREFIX = $(HOME)/.local
LINK = $(CC)
$L: $O; ar rv $@ $O && ranlib $@
.c.o:; $(CC) -c $(CFLAGS) -w -I. $<
install: $L
	mkdir -p $(PREFIX)/lib && \
	cp $L $(PREFIX)/lib

clean:; rm -f -- $L $O

include dep.mk
