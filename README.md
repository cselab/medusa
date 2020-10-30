# Install

Needs `mpifort` and `fftw3f`.

```
$ MAKEFLAGS=-j ./install.sh
```

# Run

Configuration file: [medusa/Ctrl](medusa/Ctrl). Dumps default*.nc files
file.

```
cd medusa
mpirun -n 8 ./ppm_tv
```
