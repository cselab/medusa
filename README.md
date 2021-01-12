# Install

Needs `netcdf` and `mpif90`.

```
$ MAKEFLAGS=-j make
```

Compilers and flags are in [conf.mk](conf.mk).

falcon:
```
$ module load gnu openmpi && MAKEFLAGS=-j make
````

# Run

Configuration file: [medusa/Ctrl](medusa/Ctrl). Dumps pvti files.

```
$ cd medusa
$ mpiexec -n 8 ./ppm_tv
```
