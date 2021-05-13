cmake_minimum_required(VERSION 3.14)

project(SCENE
  VERSION 2.0.0
  DESCRIPTION "Tokamak equilibrium solver"
  LANGUAGES Fortran
  HOMEPAGE_URL "https://gitlab.com/MrBehavin/SCENEv2"
  )

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

option(SCENE_ENABLE_DOUBLE "Promotes precisions of real and complex to double" ON)

add_subdirectory(lib/dierckx)

add_subdirectory(trunk)

include(CTest)
if (BUILD_TESTING)
  add_subdirectory(tests)

  add_custom_target(check
    COMMAND ctest --output-on-failure)
endif()
