if (SCENE_DOWNLOAD_GHOST)
  message(STATUS "Downloading and configuring GHOST...")
  include(FetchContent)
  FetchContent_Declare(
    ghost
    GIT_REPOSITORY https://github.com/ZedThree/GHOST
    GIT_TAG        "cba30b43bdcc73fb87cff0724337a7d3a1bd7812"
    )
  FetchContent_MakeAvailable(ghost)
  message(STATUS "Finished configuring GHOST")
  target_compile_options(ghost PRIVATE
    $<$<Fortran_COMPILER_ID:GNU>:
      -Wno-align-commons >
      )
else()
  find_package(GHOST REQUIRED)
endif()

if (SCENE_DOWNLOAD_NEASYF)
  message(STATUS "Downloading and configuring neasy-f...")
  include(FetchContent)
  FetchContent_Declare(
    neasyf
    GIT_REPOSITORY https://github.com/PlasmaFAIR/neasy-f
    GIT_TAG        main
    )
  FetchContent_MakeAvailable(neasyf)
  message(STATUS "Finished configuring neasy-f")
else()
  find_package(neasyf REQUIRED)
endif()

find_package(netCDFFortran REQUIRED)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

add_executable(scene)

if (SCENE_HAS_NBEAMS)
  set(SCENE_NBI_SOURCE_FILE nbicur.f90)
  find_package(nbeams REQUIRED)
  target_link_libraries(scene PRIVATE nbeams::nbeams)
else()
  set(SCENE_NBI_SOURCE_FILE nbicur_no_nbeams.f90)
endif()

set(SCENE_SOURCES
  astradat.f90
  baleq.f90
  balloon.f90
  balpar.f90
  bandat.f90
  bandit.f90
  caxe_out.f90
  chease_dat.f90
  coldat.f90
  elite_data.f90
  equil.f90
  error_msg.f90
  extj.f90
  fastbs.f90
  ffdgen.f90
  flxav.f90
  flxorb.f90
  flxsur.f90
  geqdsk.f90
  getdata.f90
  gfile.f90
  ghost_interface.f90
  graphs.f90
  gs2.f90
  helena.f90
  hirsh.f90
  hirsig.f90
  init.f90
  lapack_interface.f90
  lynton.f90
  mercier.f90
  outcdf.f90
  output.f90
  param.f90
  peqdsk.f90
  profs.f90
  scene.f90
  setnt.f90
  signeo.f90
  splines.f90
  test_profs.f90
  tokameq.f90
  toms790.f90
  torcur.f90
  usrcal.f90
  ${SCENE_NBI_SOURCE_FILE}
  )

target_sources(scene PRIVATE ${SCENE_SOURCES})
target_link_libraries(scene PRIVATE netCDF::netcdff neasyf::neasyf)
target_link_libraries(scene PRIVATE ghost::ghost)
target_link_libraries(scene PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(scene PRIVATE surfit)
set_target_properties(scene PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (SCENE_ENABLE_DOUBLE)
  target_compile_options(scene PUBLIC
    $<$<Fortran_COMPILER_ID:GNU>:
      -fdefault-real-8 -fdefault-double-8 >
    $<$<Fortran_COMPILER_ID:Intel>:
      -r8 >
    )
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_options(scene PRIVATE
    $<$<Fortran_COMPILER_ID:GNU>:
      -g -Wall -fimplicit-none -fbounds-check -Wimplicit-interface -Wimplicit-procedure>
    $<$<Fortran_COMPILER_ID:Intel>:
      -g -implicitnone -warn all -nogen-interfaces -CB -traceback >
    )
endif()

# Enforce 2008 standard
target_compile_options(scene PRIVATE
  $<$<Fortran_COMPILER_ID:GNU>:
    -std=f2008 >
  )

include(Sanitizers)
enable_sanitizers(scene)
