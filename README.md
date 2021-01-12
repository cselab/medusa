# Install

Needs `netcdf` and `mpif90`.

```
$ MAKEFLAGS=-j make
```

Compilers and flags are in [conf.mk](conf.mk).

# Run

Configuration file: [medusa/Ctrl](medusa/Ctrl). Dumps netcdf and vtk
files.

```
$ cd medusa
$ mpiexec -n 8 ./ppm_tv
```
