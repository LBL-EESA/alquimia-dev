# Alquimia Developer Guide for AI Agents

This repository contains the Alquimia interface library, which provides a unified interface to various chemistry engines (PFlotran, CrunchFlow) for environmental applications. It uses a CMake build system and mixes C, C++, and Fortran.

## 1. Build Instructions

The project uses CMake. To build the project:

```bash
mkdir build
cd build
cmake ..
make -j$(nproc)
```

**Key CMake Options:**
- `ALQUIMIA_SUPERBUILD`: Set to `ON` to automatically download and build dependencies (PETSc, etc.).
- `XSDK_WITH_PFLOTRAN`: Set to `ON` to enable PFlotran support.
- `XSDK_WITH_CRUNCHFLOW`: Set to `ON` to enable CrunchFlow support.
- `BUILD_SHARED_LIBS`: Defaults to `ON`.

**Dependencies:**
- PETSc (Required)
- MPI (Required)
- PFlotran (Optional)
- CrunchFlow (Optional)

## 2. Testing

Tests are managed by CTest.

**Running Tests:**
- **Run all tests:**
  ```bash
  cd build
  ctest
  ```
- **Run a specific test:**
  Use the `-R` flag with a regex pattern.
  ```bash
  ctest -R test_alquimia_c_utils_exe
  ctest -R pflotran  # Run all PFlotran-related tests
  ```
- **List all available tests:**
  ```bash
  ctest -N
  ```
- **Verbose output (for debugging failures):**
  ```bash
  ctest -V -R <test_name>
  ```

**Common Test Targets:**
- `test_alquimia_c_utils_exe`: Unit tests for C utility functions.
- `batch_chem_*`: Batch chemistry driver tests.
- `transport_*`: Transport driver tests.

## 3. Code Style & Conventions

The codebase generally follows the **Google C++ Style Guide**, but with C idioms where appropriate.

### General
- **Indentation:** 2 spaces. No tabs.
- **File Encoding:** UTF-8.
- **Line Length:** Try to keep lines under 80-100 characters, though not strictly enforced.

### C / C++
- **Naming:**
  - **Functions:** `CamelCase` (e.g., `CreateAlquimiaInterface`, `AllocateAlquimiaState`).
  - **Variables:** `snake_case` (e.g., `engine_name`, `num_primary`).
  - **Types/Structs:** `CamelCase` (e.g., `AlquimiaInterface`, `AlquimiaEngineStatus`).
  - **Constants/Enums:** `kCamelCase` (e.g., `kAlquimiaNoError`, `kAlquimiaStringPFloTran`).
  - **Macros:** `UPPER_SNAKE_CASE` (e.g., `ALQUIMIA_HAVE_PFLOTRAN`).

- **Memory Management:**
  - Explicitly use `Allocate<Type>` and `Free<Type>` functions provided by the library (e.g., `AllocateAlquimiaState`, `FreeAlquimiaState`).
  - Always pair allocations with frees to prevent leaks.

- **Error Handling:**
  - Most functions return `void` but accept a pointer to `AlquimiaEngineStatus`.
  - Check `status->error != kAlquimiaNoError` after function calls.
  - Set helpful error messages in `status->message`.

### Fortran
- **Style:** Modern Fortran (free form).
- **Naming:**
  - **Modules:** `CamelCase_module` (e.g., `AlquimiaContainers_module`).
  - **Types:** `CamelCase`.
  - **Constants:** `kCamelCase`.
- **Interoperability:** Use `iso_c_binding` for all structs and functions exposed to C/C++.

## 4. Project Structure

- `alquimia/`: Core library source code.
  - `alquimia_interface.c/h`: Main interface definition.
  - `alquimia_containers.F90`: Fortran data structure definitions.
  - `alquimia_memory.c/h`: Memory management routines.
- `drivers/`: Example drivers and executables (e.g., `BatchChemDriver`).
- `benchmarks/`: Benchmark cases and input files.
- `unit_tests/`: CTest unit tests.
- `cmake/`: CMake modules and templates.

## 5. Agent Guidelines

- **Search First:** Always use `grep` or `glob` to find relevant files before creating new ones. The project structure is specific.
- **Respect Context:** When modifying files, look at the existing code in that file. Mimic its indentation, commenting style, and variable naming conventions exactly.
- **Safety:** Do not commit secrets. Do not change `CMakeLists.txt` logic unless explicitly requested.
- **Verification:** After making changes, try to compile the project (`make`) and run relevant tests (`ctest -R <test_name>`).
