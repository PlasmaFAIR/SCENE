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
  if(NX_CONFIG_OUTPUT)
    set(${VAR} ${NX_CONFIG_OUTPUT} PARENT_SCOPE)
  endif()
endfunction()

# Copy FILENAME from source directory to build directory
macro(scene_copy_file FILENAME)
  configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
      ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}
      COPYONLY)
endmacro()


function(scene_add_test TESTNAME TESTFILE)
  set(multiValueArgs ARGUMENTS FILES)
  cmake_parse_arguments(SCENE_TEST "" "" "${multiValueArgs}" ${ARGN})

  scene_copy_file("${TESTFILE}")

  if (SCENE_TEST_FILES)
    foreach (FILE IN ITEMS ${SCENE_TEST_FILES})
      scene_copy_file("${FILE}")
    endforeach()
  endif()

  add_test(NAME ${TESTNAME}
    COMMAND ./${TESTFILE} ${SCENE_TEST_ARGUMENTS})
endfunction()
