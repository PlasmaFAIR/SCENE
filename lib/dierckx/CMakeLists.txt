add_library(surfit)

set(SURFIT_SOURCES
  bispev.f
  fpback.f
  fpbisp.f
  fpbspl.f
  fpdisc.f
  fpgivs.f
  fporde.f
  fprank.f
  fprati.f
  fprota.f
  fpsurf.f
  surfit.f
  )

target_sources(surfit PRIVATE ${SURFIT_SOURCES})

if (SCENE_ENABLE_DOUBLE)
  target_compile_options(surfit PUBLIC
    $<$<Fortran_COMPILER_ID:GNU>:
      -fdefault-real-8 -fdefault-double-8 >
    $<$<Fortran_COMPILER_ID:Intel>:
      -r8 >
    )
endif()

