# Continuous integration

The main workflow in `workflows/dev.yml` separates dependency preparation from
the two Alquimia build paths:

1. `prepare-petsc` restores or builds a prefix installation of PETSc 3.20.0.
2. `build` compiles the chemistry libraries and runs the Alquimia test suite.
3. `standalone-build` exercises the documented superbuild option that also
   builds the PFLOTRAN and CrunchFlow executables.

The two build jobs depend on `prepare-petsc` and restore their own copy of its
cache because GitHub-hosted jobs run on separate clean machines. PFLOTRAN,
CrunchFlow, and Alquimia are deliberately rebuilt in each applicable job so
source and integration changes cannot be hidden by stale build products.

## PETSc cache

The cache contains only the installed PETSc prefix at
`/home/runner/petsc-3.20.0`. Its exact key records the Ubuntu release, CPU
architecture, PETSc version, compiler family, MPI implementation, HDF5 use,
and a manual configuration revision. Partial restore keys are intentionally not
used because PETSc binaries are sensitive to all of those choices.

On a cache miss, `prepare-petsc` builds PETSc and saves it immediately. The two
dependent jobs then restore that key. On a cache hit, PETSc compilation is
skipped. System packages are still installed in every job because the cached
shared libraries depend on the runner's OpenMPI and HDF5 libraries.

Cache entries are immutable. When changing the PETSc configure command or an
ABI-relevant system dependency, increment the final revision in
`PETSC_CACHE_KEY` (for example, from `v1` to `v2`). When changing PETSc itself,
update `PETSC_VERSION`, `PETSC_PREFIX`, and `PETSC_CACHE_KEY` together.

The workflow must remain able to recreate PETSc whenever the cache is absent
or has been evicted. Do not put credentials or other sensitive data in the
cached directory.
