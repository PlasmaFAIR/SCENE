# SCENE

SCENE is a tokamak equilibrium solver which can generate equilibria in
a variety of file formats.


## Installation

SCENE requires:

- a F2008 Fortran compiler (tested with gfortran)
- netCDF Fortran wrapper (https://github.com/Unidata/netcdf-fortran)
- LAPACK and BLAS
- the GHOST Fortran graphics library
  (https://github.com/ZedThree/GHOST currently a private repo)
- CMake

To compile SCENE, you either need to tell CMake about the location of
GHOST, or you can tell CMake to download and build it automatically:

```bash
# Configure the build, telling CMake to download GHOST
cmake . -B build -DSCENE_DOWNLOAD_GHOST=ON
# OR
# Configure the build, telling CMake where to find GHOST
cmake . -B build -DGHOST_ROOT=/path/to/ghost
# Compile SCENE
cmake --build build
# Run the tests
cmake --build build --target check
```

This will create a `scene` executable under `build/bin/`.

You also need to tell CMake about the location of the other
libraries. You can do that through the below variables using the
syntax `-D<variable>=/path/to/library`:

- NetCDF: `netCDFFortran_ROOT` for the Fortran API, and `netCDF_ROOT`
  for the base C library.
    - Note: for some operating systems that put Fortran modules in an
      related location to the headers/libraries, such as Fedora, you
      may also need to set `netCDF_Fortran_INCLUDE_DIR` to the
      directory containing `netcdf.mod`
- BLAS: `BLAS_ROOT`
- LAPACK: `LAPACK_ROOT`
