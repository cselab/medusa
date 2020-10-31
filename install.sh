#!/bin/sh

set -e

(cd netcdf
 ./configure &&
     make 'FC = mpifort' 'FCFLAGS = -O2 -g -fallow-argument-mismatch' &&
     make install 'prefix = $(HOME)/.local')

(cd metis
 make install)

(cd vizing
 make install)

(cd ppm
 make install)

(cd medusa
 make)
