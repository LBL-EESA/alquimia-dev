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
option(TPL_PFLOTRAN_LIBRARIES "List of absolute paths to PFlotran link libraries [].")
option(TPL_PFLOTRAN_INCLUDE_DIRS "List of absolute paths to PFlotran include directories [].")
option(XSDK_WITH_CRUNCHFLOW "Enables support for the CrunchFlow chemistry engine [OFF]." OFF)
option(TPL_CRUNCHFLOW_LIBRARIES "List of absolute paths to CrunchFlow link libraries [].")
option(TPL_CRUNCHFLOW_INCLUDE_DIRS "List of absolute paths to CrunchFlow include directories [].")

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
if (NOT CMAKE_CXX_COMPILER)
  if (NOT $ENV{CXX} STREQUAL "")
    set(CMAKE_CXX_COMPILER $ENV{CXX})
  else()
    set(CMAKE_CXX_COMPILER c++)
  endif()
endif()
if (NOT CMAKE_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS $ENV{CXX_FLAGS})
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
enable_language(CXX)
enable_language(Fortran)

# We declare the project here.
project (alquimia)

message("-- C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})")
message("-- Fortran compiler is ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")

# Version numbers.
set (ALQUIMIA_MAJOR_VERSION 0)
set (ALQUIMIA_MINOR_VERSION 4)
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
  set(DYLIB_SUFFIX "dylib")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
else ()
  if (LINUX EQUAL 1)
    set(SYS_FLAGS "-DLINUX=1")
    set(DYLIB_SUFFIX "so")
  else()
    if (WIN32 EQUAL 1)
      set(ALQUIMIA_HAVE_BOOL 0) # MS doesn't have reasonable C compilers.
      set(SYS_FLAGS "-DWINDOWS=1")
      set(DYLIB_SUFFIX "dll")
    endif()
  endif ()
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

  message("-- Enabled support for PFlotran chemistry engine.")
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

  message("-- Enabled support for CrunchFlow chemistry engine.")
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
  if ($ENV{PETSC_ARCH} STREQUAL "")
    message(FATAL_ERROR "PETSC_ARCH must be set for the requested engines to be enabled.")
  endif()
  if (EXISTS $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/conf/PETScConfig.cmake)
    include($ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/conf/PETScConfig.cmake)
  else()
    if (EXISTS $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/petsc/conf/PETScConfig.cmake)
      include($ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/petsc/conf/PETScConfig.cmake)
    else()
      message(FATAL_ERROR "Could not find PETSc configuration info in $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}.")
    endif()
  endif()
  add_library(petsc UNKNOWN IMPORTED)
  if (BUILD_SHARED_LIBS)
    set(PETSC_LIB_NAME libpetsc${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    set(PETSC_LIB_NAME libpetsc.a)
  endif()
  if (EXISTS $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/${PETSC_LIB_NAME})
    set(PETSC_LIBRARY $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/${PETSC_LIB_NAME})
  else()
    message(FATAL_ERROR "Could not find ${PETSC_LIB_NAME} in $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}.")
  endif()
  set_target_properties(petsc PROPERTIES IMPORTED_LOCATION ${PETSC_LIBRARY})
  set(PETSC_LIBRARIES ${PETSC_LIBRARY};${PETSC_PACKAGE_LIBS})

  # PETSc 3.5.x sniffs out Valgrind if you don't ask it to, and then doesn't 
  # add the include path to its CMake-generated file. This is a bug, so we 
  # accommodate it here. A tad excrutiating.
  if (PETSC_HAVE_VALGRIND)
    file(READ $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/conf/petscvariables petscvariables)
    string(FIND ${petscvariables} "VALGRIND_INCLUDE" start_index)
    math(EXPR start_index ${start_index}+19)
    string(SUBSTRING ${petscvariables} ${start_index} -1 valgrind_include)
    string(FIND ${valgrind_include} "\n" length)
    string(SUBSTRING ${valgrind_include} 0 ${length} valgrind_include)
    if (valgrind_include) # non-blank entry
      string(REPLACE "-I" "" valgrind_include ${valgrind_include})
      message("-- Valgrind include dir is ${valgrind_include}")
      include_directories(${valgrind_include})
    endif()
  endif()
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

# Install all source headers.
install(DIRECTORY ${PROJECT_BINARY_DIR}/include/ DESTINATION include/alquimia)

# Install miscellaneous build/test files.
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/alquimia.cmake DESTINATION share/alquimia)

