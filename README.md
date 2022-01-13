# SCENE

SCENE is a tokamak equilibrium solver which can generate equilibria in
a variety of file formats.


## Installation

SCENE requires:

- a F2008 Fortran compiler (tested with gfortran)
- [neasy-f][neasyf] netCDF Fortran wrapper
    - this requires [netCDF][netcdf] and the [netCDF-Fortran][netcdff] API
- LAPACK and BLAS
- CMake

And optionally, the [GHOST Fortran graphics library][ghost].

To build SCENE, run `cmake`:

```bash
# Configure the build
cmake . -B build
# Compile SCENE
cmake --build build
```
This will create a `scene` executable under `build/bin/`.

You can run the tests like so:

```bash
cmake --build build --target check
```

You may also need to tell CMake about the location of the other
libraries. You can do that through the below variables using the
syntax `-D<variable>=/path/to/library`:

- neasy-f: `neasyf_ROOT`
    - by default, SCENE will download neasy-f automatically as required. You can
      turn this behaviour off with `-DSCENE_DOWNLOAD_NEASYF=off`
- NetCDF: `netCDFFortran_ROOT` for the Fortran API, and `netCDF_ROOT`
  for the base C library.
    - Note: for some operating systems that put Fortran modules in an
      related location to the headers/libraries, such as Fedora, you
      may also need to set `netCDF_Fortran_INCLUDE_DIR` to the
      directory containing `netcdf.mod`
- BLAS: `BLAS_ROOT`
- LAPACK: `LAPACK_ROOT`

For example, here's how to set the location of netCDF:

```bash
cmake . -B build -DnetCDFFortran_ROOT=/path/to/netcdf/fortran
```

### Using GHOST

To compile SCENE with GHOST, you can either let CMake download and
build it automatically, or give the location of an existing installation:

```bash
# Configure the build, telling CMake to download GHOST
cmake . -B build -DSCENE_USE_GHOST=ON
# OR
# Configure the build, telling CMake where to find GHOST
cmake . -B build -DSCENE_USE_GHOST=ON -DGHOST_ROOT=/path/to/ghost
```

[ghost]: https://github.com/ZedThree/GHOST
[neasyf]: https://github.com/PlasmaFAIR/neasy-f
[netcdf]: https://www.unidata.ucar.edu/software/netcdf/
[netcdff]: https://github.com/Unidata/netcdf-fortran
