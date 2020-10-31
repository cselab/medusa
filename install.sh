#!/bin/sh

set -e

(cd netcdf
 ./configure &&
     make &&
     make install 'prefix = $(HOME)/.local')

(cd metis
 make install)

(cd vizing
 make install)

(cd ppm
 make install)

(cd medusa
 make)
