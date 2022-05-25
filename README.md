# Install

Needs `netcdf` and `mpif90`. On Debian

```
sudo apt install openmpi-bin libnetcdf-mpi-dev
```

Compilers and flags are in [conf.mk](conf.mk). Buid

```
$ MAKEFLAGS=-j4 make
```

# Run

Configuration file: [medusa/Ctrl](medusa/Ctrl). Dumps pvti files.

```
$ cd medusa
$ mpiexec -n 8 ./ppm_tv
```
