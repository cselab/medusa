PREFIX = $(HOME)/.local
FC = mpif90
FCFLAGS = -O2 -g

# gfortran < 10
# FXFLAGS = -cpp

# gfortran >= 10
FXFLAGS = -cpp -fallow-argument-mismatch -Wimplicit-interface

FFTW_FCFLAGS = `pkg-config --cflags fftw3`
FFTW_LDFLAGS = `pkg-config --libs fftw3 fftw3f fftw3l`
