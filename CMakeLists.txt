cmake_minimum_required(VERSION 3.14)

project(SCENE
  VERSION 2.0.0
  DESCRIPTION "Tokamak equilibrium solver"
  LANGUAGES Fortran
  HOMEPAGE_URL "https://gitlab.com/MrBehavin/SCENEv2"
  )

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include(CMakeDependentOption)

option(SCENE_ENABLE_DOUBLE "Promotes precisions of real and complex to double" ON)
option(SCENE_USE_GHOST "Use the GHOST graphics library" OFF)
cmake_dependent_option(SCENE_DOWNLOAD_GHOST "Download and build GHOST as part of SCENE build" ON
  "SCENE_USE_GHOST" OFF)
option(SCENE_HAS_NBEAMS "Has modified version of NBEAMS with SCENE support" OFF)
option(SCENE_DOWNLOAD_NEASYF "Download and build neasy-f as part of SCENE build" ON)

add_subdirectory(lib/dierckx)

add_subdirectory(src)

include(CTest)
if (BUILD_TESTING)
  add_subdirectory(tests)

  add_custom_target(check
    COMMAND ctest --output-on-failure)
  add_dependencies(check scene)
endif()
