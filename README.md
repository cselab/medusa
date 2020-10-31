# Install

Needs `fftw3f` and `mpif90`.

```
$ MAKEFLAGS=-j ./install.sh
```

Compilers and flags are in [conf.mk](conf.mk).

# Run

Configuration file: [medusa/Ctrl](medusa/Ctrl). Dumps default*.nc files
file.

```
cd medusa
mpirun -n 8 ./ppm_tv
```
