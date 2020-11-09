PREFIX = $(HOME)/.local
FC = mpif90
FCFLAGS = -O2 -g

# gfortran < 10
# FXFLAGS = -cpp

# gfortran >= 10
FXFLAGS = -cpp -fallow-argument-mismatch -Wimplicit-interface
