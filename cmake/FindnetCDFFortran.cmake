# FindnetCDFFortran
# ----------
#
# Find the netCDF Fortran API
#
# This module uses the ``nf-config`` helper script as a hint for the
# location of the NetCDF Fortran library. It should be in your PATH.
#
# This module will define the following variables:
#
# ::
#
#   netCDFFortran_FOUND - true if netCDFFortran was found
#   netCDFFortran_VERSION - netCDFFortran version in format Major.Minor.Release
#   netCDFFortran_INCLUDE_DIRS - Location of the netCDFFortran includes
#   netCDFFortran_LIBRARIES - Required libraries
#
# This modulefFortran  You can also set the following variables:
#
# ``netCDFFortran_ROOT``
#   Specify the path to the netCDF Fortran installation to use
#
# ``netCDFFortran_DEBUG``
#   Set to TRUE to get extra debugging output

include(SCENEfunctions)

find_package(netCDF-Fortran QUIET CONFIG)
if (netCDF-Fortran_FOUND)
  set(netCDFFortran_FOUND TRUE)
  if (NOT TARGET netCDF::netcdff)
    scene_add_library_alias(netCDF::netcdff netcdff)
  endif()
  return()
endif()

find_package(netCDF REQUIRED)

find_program(NF_CONFIG "nf-config"
  PATHS "${netCDFFortran_ROOT}"
  PATH_SUFFIXES bin
  DOC "Path to netCDF Fortran config helper"
  NO_DEFAULT_PATH
  )

find_program(NF_CONFIG "nf-config"
  DOC "Path to netCDF Fortran config helper"
  )

get_filename_component(NF_CONFIG_TMP "${NF_CONFIG}" DIRECTORY)
get_filename_component(NF_CONFIG_LOCATION "${NF_CONFIG_TMP}" DIRECTORY)
if (netCDFFortran_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " NF_CONFIG_LOCATION = ${NF_CONFIG_LOCATION}")
endif()

scene_inspect_netcdf_config(NF_HINTS_INCLUDE_DIR "${NF_CONFIG}" "--includedir")
scene_inspect_netcdf_config(NF_HINTS_PREFIX "${NF_CONFIG}" "--prefix")

find_path(netCDF_Fortran_INCLUDE_DIR
  NAMES netcdf.mod
  DOC "netCDF Fortran include directories"
  HINTS
    "${netCDF_C_INCLUDE_DIR}"
    "${NF_HINTS_INCLUDE_DIR}"
    "${NF_HINTS_PREFIX}"
    "${NF_CONFIG_LOCATION}"
  PATH_SUFFIXES
    "include"
  )
if (netCDFFortran_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " netCDF_Fortran_INCLUDE_DIR = ${netCDF_Fortran_INCLUDE_DIR}"
    " NF_HINTS_INCLUDE_DIR = ${NF_HINTS_INCLUDE_DIR}"
    " NF_HINTS_PREFIX = ${NF_HINTS_PREFIX}"
    )
endif()
mark_as_advanced(netCDF_Fortran_INCLUDE_DIR)

find_library(netCDF_Fortran_LIBRARY
  NAMES netcdff
  DOC "netCDF Fortran library"
  HINTS
    "${NF_HINTS_INCLUDE_DIR}"
    "${NF_HINTS_PREFIX}"
    "${NF_CONFIG_LOCATION}"
  PATH_SUFFIXES
    "lib" "lib64"
  )
if (netCDFFortran_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " netCDF_Fortran_LIBRARY = ${netCDF_Fortran_LIBRARY}"
    " NF_HINTS_INCLUDE_DIR = ${NF_HINTS_INCLUDE_DIR}"
    " NF_HINTS_PREFIX = ${NF_HINTS_PREFIX}"
    )
endif()
mark_as_advanced(netCDF_Fortran_LIBRARY)

scene_inspect_netcdf_config(_ncxx4_version "${NF_CONFIG}" "--version")
string(REGEX REPLACE "netCDF-cxx4 \([0-9]+\)\.\([0-9]+\)\.\([0-9]+\)" "\\1" _netcdfcxx_version_major "${_ncxx4_version}")
string(REGEX REPLACE "netCDF-cxx4 \([0-9]+\)\.\([0-9]+\)\.\([0-9]+\)" "\\2" _netcdfcxx_version_minor "${_ncxx4_version}")
string(REGEX REPLACE "netCDF-cxx4 \([0-9]+\)\.\([0-9]+\)\.\([0-9]+\)" "\\3" _netcdfcxx_version_patch "${_ncxx4_version}")
set(netCDFFortran_VERSION "${_netcdf_version_major}.${_netcdf_version_minor}.${_netcdf_version_patch}${_netcdf_version_note}")
unset(_ncxx4_version)
unset(_netcdf_version_major)
unset(_netcdf_version_minor)
unset(_netcdf_version_patch)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(netCDFFortran
  REQUIRED_VARS netCDF_Fortran_LIBRARY netCDF_Fortran_INCLUDE_DIR
  VERSION_VAR netCDFFortran_VERSION)

if (netCDFFortran_FOUND)
  set(netCDFFortran_INCLUDE_DIRS "${netCDF_Fortran_INCLUDE_DIR}")
  set(netCDFFortran_LIBRARIES "${netCDF_Fortran_LIBRARY}")

  if (NOT TARGET netCDF::netcdff)
    add_library(netCDF::netcdff UNKNOWN IMPORTED)
    set_target_properties(netCDF::netcdff PROPERTIES
      IMPORTED_LINK_INTERFACE_LIBRARIES netCDF::netcdf
      IMPORTED_LOCATION "${netCDF_Fortran_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${netCDF_Fortran_INCLUDE_DIR}")
  endif ()
endif ()