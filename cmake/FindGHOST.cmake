# FindGHOST

find_library(GHOST_LIBRARIES
  NAMES ghost
  DOC "GHOST library locations"
  )

if (GHOST_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " GHOST_LIBRARIES = ${GHOST_LIBRARIES}"
    )
endif()

find_library(GRID_XGHOST_LIBRARIES
  NAMES grid_xghost
  DOC "GRID_XGHOST library locations"
  )

if (GHOST_DEBUG)
  message(STATUS "[ ${CMAKE_CURRENT_LIST_FILE}:${CMAKE_CURRENT_LIST_LINE} ] "
    " GHOST_LIBRARIES = ${GRID_XGHOST_LIBRARIES}"
    )
endif()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (GHOST DEFAULT_MSG GHOST_LIBRARIES GRID_XGHOST_LIBRARIES)

if (GHOST_FOUND AND NOT TARGET ghost::ghost)
  add_library(ghost::grid_xghost UNKNOWN IMPORTED)
  set_target_properties(ghost::grid_xghost PROPERTIES
    IMPORTED_LOCATION "${GRID_XGHOST_LIBRARIES}"
    )

  add_library(ghost::ghost UNKNOWN IMPORTED)
  set_target_properties(ghost::ghost PROPERTIES
    IMPORTED_LINK_INTERFACE_LIBRARIES ghost::grid_xghost
    IMPORTED_LOCATION "${GHOST_LIBRARIES}"
    )
endif()
