# Minimum CMake version.
cmake_minimum_required (VERSION 2.8.12)

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Options for building Alquimia. These come from the xSDK compliance rules.
option(USE_XSDK_DEFAULTS "Set to use xSDK defaults for options [ON]." ON)
option(CMAKE_INSTALL_PREFIX "Sets installation prefix [/usr/local].")
option(XSDK_ENABLE_DEBUG "Enables Debug mode builds [OFF]." OFF)
option(BUILD_SHARED_LIBS "Builds shared libraries [ON]." ON)
option(XSDK_WITH_PFLOTRAN "Enables support for the PFlotran chemistry engine [OFF]." OFF)
#option(TPL_PFLOTRAN_LIBRARIES "List of absolute paths to PFlotran link libraries [].")
#option(TPL_PFLOTRAN_INCLUDE_DIRS "List of absolute paths to PFlotran include directories [].")
option(XSDK_WITH_CRUNCHFLOW "Enables support for the CrunchFlow chemistry engine [OFF]." OFF)
#option(TPL_CRUNCHFLOW_LIBRARIES "List of absolute paths to CrunchFlow link libraries [].")
#option(TPL_CRUNCHFLOW_INCLUDE_DIRS "List of absolute paths to CrunchFlow include directories [].")

# For now, we disable shared libs on Macs.
if (APPLE)
  set(BUILD_SHARED_LIBS OFF)
endif()

if (NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX /usr/local)
endif()

# Make sure compilers are set. This must be done before enabling languages.
if (NOT CMAKE_C_COMPILER)
  if (NOT $ENV{CC} STREQUAL "")
    set(CMAKE_C_COMPILER $ENV{CC})
  else()
    set(CMAKE_C_COMPILER cc)
  endif()
endif()
if (NOT CMAKE_C_FLAGS)
  set(CMAKE_C_FLAGS $ENV{CFLAGS})
endif()
if (NOT CMAKE_Fortran_COMPILER)
  if (NOT $ENV{FC} STREQUAL "")
    set(CMAKE_Fortran_COMPILER $ENV{FC})
  else()
    set(CMAKE_Fortran_COMPILER gfortran)
  endif()
endif()
if (NOT CMAKE_Fortran_FLAGS)
  set(CMAKE_Fortran_FLAGS $ENV{FCFLAGS})
endif()
enable_language(C)
enable_language(Fortran)

# We declare the project here.
project (alquimia)

message(STATUS "C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})")
message(STATUS "Fortran compiler is ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")

if (BUILD_SHARED_LIBS)
  message(STATUS "Alquimia will be built as a shared library.")
else()
  message(STATUS "Alquimia will be built as a static library.")
endif()

# Version numbers.
set (ALQUIMIA_MAJOR_VERSION 1)
set (ALQUIMIA_MINOR_VERSION 0)
set (ALQUIMIA_PATCH_VERSION 0)
set (ALQUIMIA_VERSION "${ALQUIMIA_MAJOR_VERSION}.${ALQUIMIA_MINOR_VERSION}.${ALQUIMIA_PATCH_VERSION}")

# General C compiler flags.
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -pedantic-errors -Wextra -Werror-implicit-function-declaration")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")

  if (BUILD_SHARED_LIBS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fPIC")
  endif()

  if (LINUX EQUAL 1)
    # Counter some of GCC's more recent stinginess on Linux.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_POSIX_C_SOURCE=200809L")# -D_BSD_SOURCE")
  endif()

elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -pedantic-errors -Wextra -Werror-implicit-function-declaration -fno-builtin")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -Wno-unused-function")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

# Fortran compiler flags.
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -W -Wall -std=gnu -pedantic -Wno-unused-variable -Wno-unused-parameter")
endif()

# Figure out the system type.
set(ALQUIMIA_HAVE_BOOL 1) # All reasonable C99 compilers have this now.
if (APPLE EQUAL 1)
  set(SYS_FLAGS "-DAPPLE=1")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
else ()
  if (LINUX EQUAL 1)
    set(SYS_FLAGS "-DLINUX=1")
  else()
    if (WIN32 EQUAL 1)
      set(ALQUIMIA_HAVE_BOOL 0) # MS doesn't have reasonable C compilers.
      set(SYS_FLAGS "-DWINDOWS=1")
    endif()
  endif ()
endif ()

# Here we make sure CMake-installed binaries use the correct runpath, and 
# that the path is not stripped during installation.
if (BUILD_SHARED_LIBS)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  set(CMAKE_SKIP_INSTALL_RPATH FALSE)
else ()
  set(CMAKE_SKIP_INSTALL_RPATH TRUE)
endif ()

# Check third-party library dependencies.
set(ALQUIMIA_NEED_PETSC 0)

if (XSDK_WITH_PFLOTRAN)
  if (NOT TPL_PFLOTRAN_LIBRARIES)
    message(FATAL_ERROR "TPL_PFLOTRAN_LIBRARIES option be set for PFlotran support to be enabled.")
  endif()
  foreach(lib ${TPL_PFLOTRAN_LIBRARIES})
    if (NOT EXISTS ${lib})
      message(FATAL_ERROR "PFlotran library not found: ${lib}")
    endif()
  endforeach()

  if (NOT TPL_PFLOTRAN_INCLUDE_DIRS)
    message(FATAL_ERROR "TPL_PFLOTRAN_INCLUDE_DIRS option be set for PFlotran support to be enabled.")
  endif()
  foreach(dir ${TPL_PFLOTRAN_INCLUDE_DIRS})
    if (NOT EXISTS ${dir})
      message(FATAL_ERROR "PFlotran include directory not found: ${dir}")
    endif()
  endforeach()

  message(STATUS "Enabled support for PFlotran chemistry engine.")
  list(APPEND ALQUIMIA_ENGINES pflotran)
  set(ALQUIMIA_HAVE_PFLOTRAN 1)
  set(ALQUIMIA_NEED_PETSC 1)
else()
  set(ALQUIMIA_HAVE_PFLOTRAN 0)
endif()

if (XSDK_WITH_CRUNCHFLOW)
  if (NOT TPL_CRUNCHFLOW_LIBRARIES)
    message(FATAL_ERROR "TPL_CRUNCHFLOW_LIBRARIES option be set for CrunchFlow support to be enabled.")
  endif()
  foreach(lib ${TPL_CRUNCHFLOW_LIBRARIES})
    if (NOT EXISTS ${lib})
      message(FATAL_ERROR "CrunchFlow library not found: ${lib}")
    endif()
  endforeach()

  if (NOT TPL_CRUNCHFLOW_INCLUDE_DIRS)
    message(FATAL "TPL_CRUNCHFLOW_INCLUDE_DIRS option be set for CrunchFlow support to be enabled.")
  endif()
  foreach(dir ${TPL_CRUNCHFLOW_INCLUDE_DIRS})
    if (NOT EXISTS ${dir})
      message(FATAL "CrunchFlow include directory not found: ${dir}")
    endif()
  endforeach()

  message(STATUS "Enabled support for CrunchFlow chemistry engine.")
  list(APPEND ALQUIMIA_ENGINES crunchflow)
  set(ALQUIMIA_HAVE_CRUNCHFLOW 1)
  set(ALQUIMIA_NEED_PETSC 1)
else()
  set(ALQUIMIA_HAVE_CRUNCHFLOW 0)
endif()

# If we're not building with any engines, there's no point!
if (NOT ALQUIMIA_ENGINES)
  message(FATAL_ERROR "ERROR: At least one chemistry engine must be enabled.")
endif()

# If needed, Check for PETSc and set things up.
if (ALQUIMIA_NEED_PETSC)
  if ($ENV{PETSC_DIR} STREQUAL "")
    message(FATAL_ERROR "PETSC_DIR must be set for the requested engines to be enabled.")
  endif()
  set(PETSC_DIR $ENV{PETSC_DIR})
  if(DEFINED ENV{PETSC_ARCH})
    set(PETSC_ARCH $ENV{PETSC_ARCH})
  endif()

  message(STATUS "--- SM: PETSC_DIR >>> ${PETSC_DIR}")
  message(STATUS "--- SM: PETSC_ARCH >>> ${PETSC_ARCH}")
  
  # Do we have settings for PETSc's libraries/includes?
  if (TPL_PETSC_INCLUDE_DIRS)
    if (NOT TPL_PETSC_LDFLAGS)
      message(FATAL_ERROR "TPL_PETSC_LDFLAGS must be given if TPL_PETSC_INCLUDE_DIRS is given.")
    endif()
  elseif (TPL_PETSC_LDFLAGS)
    message(FATAL_ERROR "TPL_PETSC_INCLUDE_DIRS must be given if TPL_PETSC_LDFLAGS is given.")
  endif()

  # Find PETSc.
  if (TPL_PETSC_INCLUDE_DIRS)

    # Validate the include directories we've been given.
    foreach(inc ${TPL_PETSC_INCLUDE_DIRS})
      if (NOT EXISTS ${inc})
        message(FATAL_ERROR "Invalid PETSc include dir: ${inc}")
      endif()
    endforeach()

    message(STATUS "Using PETSc includes: ${TPL_PETSC_INCLUDE_DIRS}")
    message(STATUS "Using PETSc library flags: ${TPL_PETSC_LDFLAGS}")
    set(PETSC_INCLUDES ${TPL_PETSC_INCLUDE_DIRS})
    set(PETSC_LDFLAGS ${TPL_PETSC_LDFLAGS})
    set(PETSC_FOUND TRUE)
  else()
    # Try to detect PETSc if we haven't been told where it is.
    find_package(PETSc)
    message(STATUS "--- SM: finding PETSc with find_package")
    message(STATUS "--- SM: PETSC_INCLUDES >>> ${PETSC_INCLUDES}")
    message(STATUS "--- SM: PETSC_LDFLAGS >>> ${PETSC_LDFLAGS}")
  endif()
  if (NOT PETSC_FOUND)
    message(FATAL_ERROR "PETSc was not found.")
  endif()
  include_directories(${PETSC_INCLUDES})

  # Include PETSc in the rpath.
  if(DEFINED ENV{PETSC_ARCH})
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:${PETSC_DIR}/${PETSC_ARCH}/lib")
  else()
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:${PETSC_DIR}/lib")
  endif()
  message(STATUS "--- SM: cmake install rpath")
endif()

# Include the binary directory in the header file search path,
# since it's where we place generated files.
include_directories("${PROJECT_BINARY_DIR}")

# Source code itself.
include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(alquimia)

# Unit testing.
enable_testing()
add_subdirectory(unit_tests)

# Drivers for benchmarks.
add_subdirectory(drivers)

# Benchmarks.
add_subdirectory(benchmarks)

# Now that we have gathered all our libraries, generate an alquimia.cmake 
# file that contains all the vital information.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Templates/alquimia.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/alquimia.cmake"
  @ONLY
)

# Install miscellaneous build/test files.
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/alquimia.cmake DESTINATION share/alquimia)

