#!/bin/sh

set -e

(cd netcdf
 ./configure --disable-f03 --disable-dependency-tracking --disable-shared &&
     make &&
     make install 'prefix = $(HOME)/.local')

(cd ppm
 make install)

(cd medusa
 make)
