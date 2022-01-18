cmake_minimum_required(VERSION 3.14)

project(SCENE
  VERSION 2.0.0
  DESCRIPTION "Tokamak equilibrium solver"
  LANGUAGES Fortran
  HOMEPAGE_URL "https://gitlab.com/MrBehavin/SCENEv2"
  )

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

##################################################
# Build config options

include(CMakeDependentOption)

option(SCENE_ENABLE_DOUBLE "Promotes precisions of real and complex to double" ON)
option(SCENE_USE_GHOST "Use the GHOST graphics library" OFF)
cmake_dependent_option(SCENE_DOWNLOAD_GHOST "Download and build GHOST as part of SCENE build" ON
  "SCENE_USE_GHOST;NOT GHOST_DIR" OFF)
option(SCENE_USE_NBEAMS "Has modified version of NBEAMS with SCENE support" OFF)
option(SCENE_DOWNLOAD_NEASYF "Download and build neasy-f as part of SCENE build" ON)
option(SCENE_DOWNLOAD_FUUID4 "Download and build fuuid4 as part of SCENE build" ON)
option(SCENE_UPDATE_GIT_SUBMODULE "Check submodules are up-to-date during build" ON)

##################################################
# Targets

include(UpdateSubmodules)
scene_update_submodules()

add_subdirectory(externals/fortran-git-version)
add_subdirectory(lib/dierckx)

add_subdirectory(src)

include(CTest)
if (BUILD_TESTING)
  add_subdirectory(tests)

  add_custom_target(check
    COMMAND ctest --output-on-failure)
  add_dependencies(check scene)
endif()
