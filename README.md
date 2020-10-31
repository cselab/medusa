# Install

Needs `mpifort` and `fftw3f`.

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
