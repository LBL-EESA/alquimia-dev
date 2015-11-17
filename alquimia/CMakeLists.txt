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
if (XSDK_WITH_PFLOTRAN)
  include_directories(${TPL_PFLOTRAN_INCLUDE_DIRS})
  list(APPEND alquimia_sources pflotran_alquimia_interface.F90 pflotran_alquimia_wrappers.F90)
endif()

# CrunchFlow engine.
if (XSDK_WITH_CRUNCHFLOW)
  include_directories(${TPL_CRUNCHFLOW_INCLUDE_DIRS})
  list(APPEND alquimia_sources crunch_alquimia_interface.F90 crunch_alquimia_wrappers.F90)
endif()

add_library(alquimia ${ALQUIMIA_SOURCES})

# Install library.
install(TARGETS alquimia DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/alquimia.h)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/alquimia_version.h)
install(FILES ${HEADERS} DESTINATION include/alquimia)
