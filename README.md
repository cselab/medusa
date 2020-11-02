# Install

Needs `fftw3f` and `mpif90`.

```
$ MAKEFLAGS=-j ./install.sh
```

Compilers and flags are in [conf.mk](conf.mk).

# Run

Configuration file: [medusa/Ctrl](medusa/Ctrl). Dumps netcdf and vtk
files.

```
cd medusa
mpirun -n 8 ./ppm_tv
```
