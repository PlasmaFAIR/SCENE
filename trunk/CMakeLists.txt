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
else()
  find_package(GHOST REQUIRED)
endif()

find_package(netCDFFortran REQUIRED)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

add_executable(scene)

set(SCENE_SOURCES
  astradat.f90
  baleq.f90
  balloon.f90
  balpar.f90
  bandat.f90
  bandit.f90
  beamconsts.f
  beams_mfp.f
  calcbeams.f
  caxe_out.f90
  chease_dat.f90
  coldat.f90
  coulomb.f
  eisplit.f
  elite_data.f90
  equil.f90
  erf.f
  error_msg.f90
  extj.f90
  fastbs.f90
  fastions.f
  ffdgen.f90
  flxav.f90
  flxorb.f90
  flxsur.f90
  frate.f
  gausswts.f
  geqdsk.f90
  getdata.f90
  getrho.f
  gfile.f90
  ghost_interface.f90
  graphs.f90
  gs2.f90
  helena.f90
  hirsh.f90
  hirsig.f90
  hofr.f
  hunt.f
  init.f90
  lynton.f90
  mercier.f90
  nbicur.f90
  outcdf.f90
  output.f90
  param.f90
  peqdsk.f90
  profs.f90
  qsimp.f
  rinteg.f
  scene.f90
  setnt.f90
  sig_olson.f
  sigfit.f
  signeo.f90
  sinteg.f
  splines.f
  ssum.f
  test_profs.f90
  tokameq.f90
  toms790.f90
  torcur.f90
  usrcal.f90
  zinteg.f
  )

target_sources(scene PRIVATE ${SCENE_SOURCES})
target_link_libraries(scene ghost::ghost netCDF::netcdff)
target_link_libraries(scene ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(scene surfit)
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
      -g -Wall -fimplicit-none -fbounds-check >
    $<$<Fortran_COMPILER_ID:Intel>:
      -g -implicitnone -warn all -nogen-interfaces -CB -traceback >
    )
endif()

# Enforce 2008 standard
target_compile_options(scene PRIVATE
  $<$<Fortran_COMPILER_ID:GNU>:
    -std=f2008 >
  )
