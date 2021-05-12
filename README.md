# SCENE

SCENE is a tokamak equilibrium solver which can generate equilibria in
a variety of file formats.


## Installation

SCENE requires:

- a F2008 Fortran compiler (tested with gfortran)
- netCDF Fortran wrapper (https://github.com/Unidata/netcdf-fortran)
- the GHOST Fortran graphics library
  (https://github.com/ZedThree/GHOST currently a private repo)

To compile SCENE, you need to set two variables that point at the
install locations of the above two libraries:

```bash
  make YORK_GHOST_DIR=/path/to/ghost/install \
       YORK_NETCDF_DIR=/path/to/netcdf/install
```

This will create a `scene` executable under `trunk/`
