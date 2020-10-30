.POSIX:
.SUFFIXES: .o .cpp
L = libvizing.a
all: $L
O = vizing_coloring.o

PREFIX = $(HOME)/.local
CXXFLAGS = -Ofast -g
CXX = g++
LINK = $(CXX)
$L: $O; ar rv $@ $O && ranlib $@
.cpp.o:; $(CXX) -c $(CXXFLAGS) $<
install: $L
	mkdir -p $(PREFIX)/lib && \
	cp $L $(PREFIX)/lib

clean:; rm -f -- $L $O

vizing_coloring: \
adj.h\
Array.h\
Boolean.h\
color.h\
edge.h\
graph.h\
list.h\
