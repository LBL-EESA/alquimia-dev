# Build a configuration header file from our options.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/alquimia.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/alquimia.h"
)

# Core library sources.
set(ALQUIMIA_SOURCES alquimia.c;alquimia_constants.c;alquimia_interface.c;
    alquimia_memory.c;alquimia_util.c;alquimia_containers.F90;
    c_f_interface_module.F90;alquimia_fortran_interface_mod.F90)

# PFlotran engine.
if (ALQUIMIA_HAVE_PFLOTRAN)
  list(APPEND ALQUIMIA_TPLS ${TPL_PFLOTRAN_LIBRARIES})
  list(APPEND ALQUIMIA_INCLUDE_DIRS ${TPL_PFLOTRAN_INCLUDE_DIRS})
  list(APPEND ALQUIMIA_SOURCES pflotran_alquimia_interface.F90 pflotran_alquimia_wrappers.F90)
  add_library(pflotran_petsc_stubs pflotran_petsc_stubs.F90)
  list(APPEND ALQUIMIA_TPLS pflotran_petsc_stubs)
endif()

# CrunchFlow engine.
if (ALQUIMIA_HAVE_CRUNCHFLOW)
  list(APPEND ALQUIMIA_TPLS ${TPL_CRUNCHFLOW_LIBRARIES})
  list(APPEND ALQUIMIA_INCLUDE_DIRS ${TPL_CRUNCHFLOW_INCLUDE_DIRS})
  list(APPEND ALQUIMIA_SOURCES crunch_alquimia_interface.F90 crunch_alquimia_wrappers.F90)
endif()

if (ALQUIMIA_NEED_PETSC)
  list(APPEND ALQUIMIA_INCLUDE_DIRS $ENV{PETSC_DIR}/include)
  list(APPEND ALQUIMIA_INCLUDE_DIRS $ENV{PETSC_DIR}/include/petsc)
  list(APPEND ALQUIMIA_INCLUDE_DIRS $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include)
  if (PETSC_LIBRARIES)
    list(APPEND ALQUIMIA_TPLS ${PETSC_LIBRARIES})
  endif()
  if (PETSC_LDFLAGS)
    list(APPEND ALQUIMIA_TPLS ${PETSC_LDFLAGS})
  endif()
endif()

include_directories(${ALQUIMIA_INCLUDE_DIRS})
add_library(alquimia ${ALQUIMIA_SOURCES})
target_link_libraries(alquimia ${ALQUIMIA_TPLS})

# Stuff for provenance.
if (BUILD_SHARED_LIBS)
  set(ALQUIMIA_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/libalquimia${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(ALQUIMIA_LIBRARY ${CMAKE_INSTALL_PREFIX}/lib/libalquimia.a)
endif()
set(ALQUIMIA_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include/alquimia;${ALQUIMIA_INCLUDE_DIRS} PARENT_SCOPE)
set(ALQUIMIA_TPLS ${ALQUIMIA_TPLS} PARENT_SCOPE)
set(ALQUIMIA_LIBRARIES ${ALQUIMIA_LIBRARY};${ALQUIMIA_TPLS} PARENT_SCOPE)

# Install library.
if (ALQUIMIA_HAVE_PFLOTRAN)
  install(TARGETS alquimia pflotran_petsc_stubs DESTINATION lib)
else()
  install(TARGETS alquimia DESTINATION lib)
endif()

# Install headers, including generated files.
file(GLOB HEADERS *.h)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/alquimia.h)
install(FILES ${HEADERS} DESTINATION include/alquimia)
