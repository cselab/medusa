PREFIX = $(HOME)/.local
FC = mpif90
FCFLAGS = -O2 -g

# gfortran < 10
#FXFLAGS = -cpp

# gfortran >= 10
FXFLAGS = -cpp -fallow-argument-mismatch

NETCDF_CFLAGS = `pkg-config --cflags netcdf-mpi`
NETCDF_LDFLAGS = `pkg-config --libs netcdf-mpi`

#NETCDF_CFLAGS = -I/usr/lib/x86_64-linux-gnu/netcdf/mpi/include
#NETCDF_LDFLAGS = -L/usr/lib/x86_64-linux-gnu/netcdf/mpi -lnetcdf
