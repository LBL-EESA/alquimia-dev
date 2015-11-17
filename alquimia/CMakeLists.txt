# Build a configuration header file from our options.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/alquimia.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/alquimia.h"
)

# Core library sources.
set(ALQUIMIA_SOURCES alquimia.c;alquimia_constants.c;alquimia_interface.c;
    alquimia_memory.c;alquimia_util.c;alquimia_containers.F90;
    c_f_interface_module.F90)

# PFlotran engine.
if (ALQUIMIA_HAVE_PFLOTRAN)
  include_directories($ENV{PETSC_DIR}/include)
  include_directories($ENV{PETSC_DIR}/include/petsc)
  include_directories($ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include)
  include_directories(${TPL_PFLOTRAN_INCLUDE_DIRS})
  list(APPEND ALQUIMIA_SOURCES pflotran_alquimia_interface.F90 pflotran_alquimia_wrappers.F90)
endif()

# CrunchFlow engine.
if (ALQUIMIA_HAVE_CRUNCHFLOW)
  include_directories(${TPL_CRUNCHFLOW_INCLUDE_DIRS})
  list(APPEND ALQUIMIA_SOURCES crunch_alquimia_interface.F90 crunch_alquimia_wrappers.F90)
endif()

link_libraries(${TPL_PFLOTRAN_LIBRARIES};${PETSC_LIBRARIES})
add_library(alquimia ${ALQUIMIA_SOURCES})

# Install library.
install(TARGETS alquimia DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/alquimia.h)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/alquimia_version.h)
install(FILES ${HEADERS} DESTINATION include/alquimia)
