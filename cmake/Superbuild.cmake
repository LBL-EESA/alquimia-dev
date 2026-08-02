include(ExternalProject)

# Set the install directory for dependencies
set(INSTALL_DIR "${CMAKE_BINARY_DIR}/install")

# Chemistry engine options
option(XSDK_WITH_PFLOTRAN "Enables support for the PFlotran chemistry engine [ON]." ON)
option(XSDK_WITH_CRUNCHFLOW "Enables support for the CrunchFlow chemistry engine [ON]." ON)
option(ALQUIMIA_BUILD_STANDALONE_ENGINES "Build standalone versions of requested chemistry engines [OFF]." OFF)

if (NOT XSDK_WITH_PFLOTRAN AND NOT XSDK_WITH_CRUNCHFLOW)
  message(FATAL_ERROR "At least one chemistry engine must be enabled (XSDK_WITH_PFLOTRAN or XSDK_WITH_CRUNCHFLOW).")
endif()

# Pass down compiler/flags to external projects
set(COMMON_CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
)

# HDF5 is needed only by standalone PFLOTRAN. When PETSc is built by this
# superbuild, let PETSc build a compatible parallel HDF5 with Fortran bindings.
set(ALQUIMIA_BUILD_STANDALONE_PFLOTRAN OFF)
set(PETSC_HDF5_ARGS --with-hdf5=0)
if(ALQUIMIA_BUILD_STANDALONE_ENGINES AND XSDK_WITH_PFLOTRAN)
  set(ALQUIMIA_BUILD_STANDALONE_PFLOTRAN ON)
  set(PETSC_HDF5_ARGS
      --download-hdf5=1
      --with-hdf5-fortran-bindings=1)
endif()

find_package(BLAS QUIET)
find_package(LAPACK QUIET)
if(BLAS_FOUND AND LAPACK_FOUND AND FALSE)
  message(STATUS "Found system BLAS/LAPACK")
  set(PETSC_BLASLAPACK_ARGS "")
else()
  message(STATUS "BLAS/LAPACK not found, will be downloaded by PETSc")
  set(PETSC_BLASLAPACK_ARGS "--download-fblaslapack=1")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(PETSC_DEBUG_ARG "--with-debugging=1")
else()
  set(PETSC_DEBUG_ARG "--with-debugging=0")
endif()

function(configure_external_petsc_hdf5 petsc_variables)
  file(STRINGS "${petsc_variables}" hdf5_variable_lines
       REGEX "^HDF5_(LIB|INCLUDE) *=")
  foreach(line ${hdf5_variable_lines})
    if(line MATCHES "^HDF5_INCLUDE *=(.*)")
      set(hdf5_include_flags "${CMAKE_MATCH_1}")
    elseif(line MATCHES "^HDF5_LIB *=(.*)")
      set(hdf5_link_flags "${CMAKE_MATCH_1}")
    endif()
  endforeach()

  string(REGEX MATCHALL "-I[^ ]+" hdf5_include_options
         "${hdf5_include_flags}")
  set(hdf5_include_dir "")
  foreach(option ${hdf5_include_options})
    string(SUBSTRING "${option}" 2 -1 candidate_include_dir)
    if(EXISTS "${candidate_include_dir}/H5pubconf.h")
      set(hdf5_include_dir "${candidate_include_dir}")
    endif()
  endforeach()

  if(hdf5_include_dir)
    file(STRINGS "${hdf5_include_dir}/H5pubconf.h" hdf5_parallel_lines
         REGEX "^ *#define +H5_HAVE_PARALLEL +1")
  endif()

  string(REGEX MATCHALL "-L[^ ]+" hdf5_library_options
         "${hdf5_link_flags}")
  set(hdf5_library_dirs "")
  foreach(option ${hdf5_library_options})
    string(SUBSTRING "${option}" 2 -1 candidate_library_dir)
    list(APPEND hdf5_library_dirs "${candidate_library_dir}")
  endforeach()

  unset(ALQUIMIA_HDF5_LIBRARY CACHE)
  unset(ALQUIMIA_HDF5_HL_LIBRARY CACHE)
  unset(ALQUIMIA_HDF5_FORTRAN_LIBRARY CACHE)
  unset(ALQUIMIA_HDF5_FORTRAN_HL_LIBRARY CACHE)
  find_library(ALQUIMIA_HDF5_LIBRARY NAMES hdf5
               PATHS ${hdf5_library_dirs} NO_DEFAULT_PATH)
  find_library(ALQUIMIA_HDF5_HL_LIBRARY NAMES hdf5_hl
               PATHS ${hdf5_library_dirs} NO_DEFAULT_PATH)
  find_library(ALQUIMIA_HDF5_FORTRAN_LIBRARY NAMES hdf5_fortran
               PATHS ${hdf5_library_dirs} NO_DEFAULT_PATH)
  find_library(ALQUIMIA_HDF5_FORTRAN_HL_LIBRARY
               NAMES hdf5hl_fortran hdf5_hl_fortran
               PATHS ${hdf5_library_dirs} NO_DEFAULT_PATH)

  if(NOT hdf5_include_dir OR NOT hdf5_parallel_lines OR
     NOT ALQUIMIA_HDF5_LIBRARY OR NOT ALQUIMIA_HDF5_HL_LIBRARY OR
     NOT ALQUIMIA_HDF5_FORTRAN_LIBRARY OR
     NOT ALQUIMIA_HDF5_FORTRAN_HL_LIBRARY)
    message(FATAL_ERROR
      "Standalone PFLOTRAN requires the provided PETSc installation to be "
      "built with parallel HDF5, including its Fortran and high-level "
      "interfaces. The HDF5 recorded in ${petsc_variables} is not suitable.")
  endif()

  get_filename_component(hdf5_library_dir
                         "${ALQUIMIA_HDF5_LIBRARY}" DIRECTORY)
  set(ALQUIMIA_HDF5_INCLUDE_DIR "${hdf5_include_dir}" PARENT_SCOPE)
  set(ALQUIMIA_HDF5_LIBRARY_DIR "${hdf5_library_dir}" PARENT_SCOPE)
  set(ALQUIMIA_HDF5_LINK_FLAGS
      "${ALQUIMIA_HDF5_FORTRAN_HL_LIBRARY} ${ALQUIMIA_HDF5_FORTRAN_LIBRARY} ${hdf5_link_flags} -lz"
      PARENT_SCOPE)
endfunction()

# Resolve PETSc from explicit CMake variables first, then the environment.
set(ALQUIMIA_PETSC_DIR "${PETSC_DIR}")
if(NOT ALQUIMIA_PETSC_DIR AND NOT "$ENV{PETSC_DIR}" STREQUAL "")
  set(ALQUIMIA_PETSC_DIR "$ENV{PETSC_DIR}")
endif()
set(ALQUIMIA_PETSC_ARCH "${PETSC_ARCH}")
if(NOT DEFINED PETSC_ARCH AND DEFINED ENV{PETSC_ARCH})
  set(ALQUIMIA_PETSC_ARCH "$ENV{PETSC_ARCH}")
endif()

if(ALQUIMIA_PETSC_DIR)
  get_filename_component(ALQUIMIA_PETSC_DIR "${ALQUIMIA_PETSC_DIR}" ABSOLUTE)
  set(PETSC_VERSION_HEADER "${ALQUIMIA_PETSC_DIR}/include/petscversion.h")
  if(NOT EXISTS "${PETSC_VERSION_HEADER}")
    message(FATAL_ERROR
      "PETSC_DIR does not contain include/petscversion.h: ${ALQUIMIA_PETSC_DIR}")
  endif()

  file(STRINGS "${PETSC_VERSION_HEADER}" PETSC_VERSION_LINES
       REGEX "^#define PETSC_VERSION_(MAJOR|MINOR|SUBMINOR) ")
  foreach(line ${PETSC_VERSION_LINES})
    if(line MATCHES "PETSC_VERSION_MAJOR +([0-9]+)")
      set(PETSC_VERSION_MAJOR "${CMAKE_MATCH_1}")
    elseif(line MATCHES "PETSC_VERSION_MINOR +([0-9]+)")
      set(PETSC_VERSION_MINOR "${CMAKE_MATCH_1}")
    elseif(line MATCHES "PETSC_VERSION_SUBMINOR +([0-9]+)")
      set(PETSC_VERSION_SUBMINOR "${CMAKE_MATCH_1}")
    endif()
  endforeach()
  if(NOT DEFINED PETSC_VERSION_MAJOR OR
     NOT DEFINED PETSC_VERSION_MINOR OR
     NOT DEFINED PETSC_VERSION_SUBMINOR)
    message(FATAL_ERROR "Could not read the PETSc version from ${PETSC_VERSION_HEADER}")
  endif()
  set(ALQUIMIA_PETSC_VERSION
      "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}")
  if(ALQUIMIA_PETSC_VERSION VERSION_LESS "3.20.0")
    message(FATAL_ERROR
      "PETSc 3.20.0 or newer is required; found ${ALQUIMIA_PETSC_VERSION} "
      "in ${ALQUIMIA_PETSC_DIR}")
  endif()

  if(ALQUIMIA_PETSC_ARCH)
    set(PETSC_VARIABLES
        "${ALQUIMIA_PETSC_DIR}/${ALQUIMIA_PETSC_ARCH}/lib/petsc/conf/petscvariables")
    set(ALQUIMIA_PETSC_LIB_DIR
        "${ALQUIMIA_PETSC_DIR}/${ALQUIMIA_PETSC_ARCH}/lib")
  else()
    set(PETSC_VARIABLES
        "${ALQUIMIA_PETSC_DIR}/lib/petsc/conf/petscvariables")
    set(ALQUIMIA_PETSC_LIB_DIR "${ALQUIMIA_PETSC_DIR}/lib")
  endif()
  if(NOT EXISTS "${PETSC_VARIABLES}")
    message(FATAL_ERROR
      "PETSC_DIR/PETSC_ARCH do not identify a usable PETSc installation: "
      "${PETSC_VARIABLES} was not found")
  endif()

  message(STATUS
    "Using external PETSc ${ALQUIMIA_PETSC_VERSION}: ${ALQUIMIA_PETSC_DIR}")
  if(ALQUIMIA_BUILD_STANDALONE_PFLOTRAN)
    configure_external_petsc_hdf5("${PETSC_VARIABLES}")
    message(STATUS
      "Using parallel HDF5 from the provided PETSc: ${ALQUIMIA_HDF5_LIBRARY_DIR}")
  endif()
  add_custom_target(petsc)
else()
  set(ALQUIMIA_PETSC_DIR "${INSTALL_DIR}")
  set(ALQUIMIA_PETSC_ARCH "")
  set(ALQUIMIA_PETSC_LIB_DIR "${INSTALL_DIR}/lib")
  if(NOT ALQUIMIA_PYTHON_EXECUTABLE)
    find_program(ALQUIMIA_PYTHON_EXECUTABLE
                 NAMES python3.12 python3.11 python3.10 python3.9 python3.8
                       python3.7 python3.6 python3.5 python3.4 python3
                 DOC "Python interpreter used to configure PETSc")
  endif()
  if(NOT ALQUIMIA_PYTHON_EXECUTABLE)
    message(FATAL_ERROR
      "A Python interpreter compatible with PETSc 3.20 was not found. Set "
      "ALQUIMIA_PYTHON_EXECUTABLE to a Python 3 interpreter with xdrlib.")
  endif()
  execute_process(
      COMMAND "${ALQUIMIA_PYTHON_EXECUTABLE}" -c
              "import sys, xdrlib; assert sys.version_info >= (3, 4)"
      RESULT_VARIABLE ALQUIMIA_PYTHON_CHECK_RESULT
      OUTPUT_QUIET
      ERROR_QUIET)
  if(NOT ALQUIMIA_PYTHON_CHECK_RESULT EQUAL 0)
    message(FATAL_ERROR
      "${ALQUIMIA_PYTHON_EXECUTABLE} cannot configure PETSc 3.20 because "
      "it does not provide Python's xdrlib module. Set "
      "ALQUIMIA_PYTHON_EXECUTABLE to a compatible Python, such as "
      "/usr/bin/python3.")
  endif()
  message(STATUS
    "Using Python to configure PETSc: ${ALQUIMIA_PYTHON_EXECUTABLE}")
  if(ALQUIMIA_BUILD_STANDALONE_PFLOTRAN)
    set(ALQUIMIA_HDF5_INCLUDE_DIR "${INSTALL_DIR}/include")
    set(ALQUIMIA_HDF5_LIBRARY_DIR "${INSTALL_DIR}/lib")
    set(ALQUIMIA_HDF5_LINK_FLAGS
        "-L${INSTALL_DIR}/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz")
  endif()
  ExternalProject_Add(petsc
      GIT_REPOSITORY https://gitlab.com/petsc/petsc.git
      GIT_TAG v3.20.0
      PREFIX ${CMAKE_BINARY_DIR}/external/petsc
      BUILD_IN_SOURCE 1
      UPDATE_DISCONNECTED 1
      CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env -- ${ALQUIMIA_PYTHON_EXECUTABLE} ./configure --prefix=${INSTALL_DIR} --with-mpi=1 ${PETSC_DEBUG_ARG} --with-shared-libraries=1 ${PETSC_HDF5_ARGS} ${PETSC_BLASLAPACK_ARGS}
      BUILD_COMMAND make
      INSTALL_COMMAND make install
  )
endif()

set(ALQUIMIA_DEPS petsc)
set(ALQUIMIA_EXTRA_ARGS)

# PFLOTRAN
if (XSDK_WITH_PFLOTRAN)
  ExternalProject_Add(pflotran
      DEPENDS petsc
      GIT_REPOSITORY https://bitbucket.org/pflotran/pflotran
      GIT_TAG v5.0.0
      PREFIX ${CMAKE_BINARY_DIR}/external/pflotran
      CONFIGURE_COMMAND ""
      UPDATE_DISCONNECTED 1
      BUILD_COMMAND make -C src/pflotran libpflotranchem.a PETSC_DIR=${ALQUIMIA_PETSC_DIR} "PETSC_ARCH=${ALQUIMIA_PETSC_ARCH}"
      BUILD_IN_SOURCE 1
      INSTALL_COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/src/pflotran/libpflotranchem.a ${INSTALL_DIR}/lib/libpflotranchem.a
              COMMAND ${CMAKE_COMMAND} -E make_directory ${INSTALL_DIR}/include/pflotran
              COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/src/pflotran ${INSTALL_DIR}/include/pflotran
  )
  list(APPEND ALQUIMIA_DEPS pflotran)
  list(APPEND ALQUIMIA_EXTRA_ARGS 
       -DXSDK_WITH_PFLOTRAN=ON
       -DTPL_PFLOTRAN_LIBRARIES=${INSTALL_DIR}/lib/libpflotranchem.a
       -DTPL_PFLOTRAN_INCLUDE_DIRS=${INSTALL_DIR}/include/pflotran)

  if (ALQUIMIA_BUILD_STANDALONE_ENGINES)
    set(PFLOTRAN_STANDALONE_EXTRA_MAKE_ARGS
        "have_hdf5=1"
        "HDF5_LIB=${ALQUIMIA_HDF5_LIBRARY_DIR}"
        "HDF5_INCLUDE=${ALQUIMIA_HDF5_INCLUDE_DIR}"
        "LIBS=${ALQUIMIA_HDF5_LINK_FLAGS}")

    ExternalProject_Add(pflotran_standalone
        DEPENDS petsc
        GIT_REPOSITORY https://bitbucket.org/pflotran/pflotran
        GIT_TAG v5.0.0
        PREFIX ${CMAKE_BINARY_DIR}/external/pflotran_standalone
        CONFIGURE_COMMAND ""
        UPDATE_DISCONNECTED 1
        BUILD_COMMAND make -C src/pflotran pflotran PETSC_DIR=${ALQUIMIA_PETSC_DIR} "PETSC_ARCH=${ALQUIMIA_PETSC_ARCH}" ${PFLOTRAN_STANDALONE_EXTRA_MAKE_ARGS}
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${INSTALL_DIR}/bin
                COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/src/pflotran/pflotran ${INSTALL_DIR}/bin/pflotran
    )
  endif()
else()
  list(APPEND ALQUIMIA_EXTRA_ARGS -DXSDK_WITH_PFLOTRAN=OFF)
endif()

# CrunchFlow
if (XSDK_WITH_CRUNCHFLOW)
  ExternalProject_Add(crunchflow
      DEPENDS petsc
      GIT_REPOSITORY https://bitbucket.org/crunchflow/crunchtope-dev
      GIT_TAG master
      PREFIX ${CMAKE_BINARY_DIR}/external/crunchflow
      CONFIGURE_COMMAND ""
      UPDATE_DISCONNECTED 1
      PATCH_COMMAND git apply --check source/MakefileForAlquimia.patch && git apply source/MakefileForAlquimia.patch || echo "Patch already applied or failed"
      BUILD_COMMAND make -C source libcrunchchem.a PETSC_DIR=${ALQUIMIA_PETSC_DIR} "PETSC_ARCH=${ALQUIMIA_PETSC_ARCH}"
      BUILD_IN_SOURCE 1
      INSTALL_COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/source/libcrunchchem.a ${INSTALL_DIR}/lib/libcrunchchem.a
              COMMAND ${CMAKE_COMMAND} -E make_directory ${INSTALL_DIR}/include/crunchflow
              COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/source ${INSTALL_DIR}/include/crunchflow
  )
  list(APPEND ALQUIMIA_DEPS crunchflow)
  list(APPEND ALQUIMIA_EXTRA_ARGS 
       -DXSDK_WITH_CRUNCHFLOW=ON
       -DTPL_CRUNCHFLOW_LIBRARIES=${INSTALL_DIR}/lib/libcrunchchem.a
       -DTPL_CRUNCHFLOW_INCLUDE_DIRS=${INSTALL_DIR}/include/crunchflow)

  if (ALQUIMIA_BUILD_STANDALONE_ENGINES)
    ExternalProject_Add(crunchflow_standalone
        DEPENDS petsc
        GIT_REPOSITORY https://bitbucket.org/crunchflow/crunchtope-dev
        GIT_TAG master
        PREFIX ${CMAKE_BINARY_DIR}/external/crunchflow_standalone
        CONFIGURE_COMMAND ""
        UPDATE_DISCONNECTED 1
        PATCH_COMMAND sed -i "s/chkopts//g" source/Makefile
        BUILD_COMMAND make -C source CrunchMain PETSC_DIR=${ALQUIMIA_PETSC_DIR} "PETSC_ARCH=${ALQUIMIA_PETSC_ARCH}"
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${INSTALL_DIR}/bin
                COMMAND ${CMAKE_COMMAND} -E copy <SOURCE_DIR>/source/CrunchTope ${INSTALL_DIR}/bin/crunchflow
    )
  endif()
else()
  list(APPEND ALQUIMIA_EXTRA_ARGS -DXSDK_WITH_CRUNCHFLOW=OFF)
endif()

# Alquimia itself
ExternalProject_Add(alquimia_core
    DEPENDS ${ALQUIMIA_DEPS}
    SOURCE_DIR ${CMAKE_SOURCE_DIR}
    BINARY_DIR ${CMAKE_BINARY_DIR}/alquimia-build
    INSTALL_DIR ${INSTALL_DIR}
    CMAKE_ARGS
        ${COMMON_CMAKE_ARGS}
        ${ALQUIMIA_EXTRA_ARGS}
        -DPETSC_DIR=${ALQUIMIA_PETSC_DIR}
        -DPETSC_ARCH=${ALQUIMIA_PETSC_ARCH}
        -DALQUIMIA_SUPERBUILD=OFF
)

# Forward the test target to the inner build
# We remove the dependency on alquimia_core so that 'make test' doesn't 
# trigger a re-check of all dependencies.
add_custom_target(test
    COMMAND ${CMAKE_COMMAND} -E env LD_LIBRARY_PATH=${INSTALL_DIR}/lib:${ALQUIMIA_PETSC_LIB_DIR}:$ENV{LD_LIBRARY_PATH} ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}/alquimia-build --target test
)

# Rule to install the contents of the local install directory to the final destination
install(DIRECTORY ${INSTALL_DIR}/ DESTINATION .)
