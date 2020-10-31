#!/bin/sh

set -e

(cd netcdf
 make install)

(cd metis
 make install)

(cd vizing
 make install)

(cd ppm
 make install)

(cd medusa
 make)
