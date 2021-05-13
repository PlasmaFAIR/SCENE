# Add an alias for an imported target
# Workaround for CMAke < 3.18
# Taken from https://github.com/conan-io/conan/issues/2125#issuecomment-351176653
function(scene_add_library_alias dst src)
  add_library(${dst} INTERFACE IMPORTED)
  foreach(name INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES INTERFACE_COMPILE_DEFINITIONS INTERFACE_COMPILE_OPTIONS)
    get_property(value TARGET ${src} PROPERTY ${name} )
    set_property(TARGET ${dst} PROPERTY ${name} ${value})
  endforeach()
endfunction()


# Call nx-config with an argument, and append the resulting path to a list
# Taken from https://github.com/LiamBindle/geos-chem/blob/feature/CMake/CMakeScripts/FindNetCDF.cmake
function(scene_inspect_netcdf_config VAR NX_CONFIG ARG)
  execute_process(
    COMMAND ${NX_CONFIG} ${ARG}
    OUTPUT_VARIABLE NX_CONFIG_OUTPUT
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(EXISTS "${NX_CONFIG_OUTPUT}")
    set(${VAR} ${NX_CONFIG_OUTPUT} PARENT_SCOPE)
  endif()
endfunction()
