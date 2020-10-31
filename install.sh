#!/bin/sh

set -e

(cd netcdf &&
./configure &&
make 'FC = mpifort' 'FFLAGS = -O2 -g -fallow-argument-mismatch' &&
make install 'prefix = $(HOME)/.local')

(cd metis
 make install)

(cd vizing
 make install)

(cd ppm
 make install)
# falcon
# (cd ppm ; make install 'FCFLAGS = -I$(FFTW_INCLUDE) -O3 -g')

(cd medusa
 make)
