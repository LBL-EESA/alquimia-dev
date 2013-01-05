#!/usr/bin/env python

from __future__ import print_function

import os

import waflib

top = '.'
build = 'build'

def unique_list(seq):
   seen = {}
   result = []
   for item in seq:
       marker = item
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

def set_fortran_stdlib(ctx):
    fortran_libs = unique_list(ctx.check_fortran_clib())
    flibs = []
    fpaths = []
    for l in fortran_libs:
        if l[0:2] == "-L":
            fpaths.append(l[2:])
        if l[0:2] == "-l":
            flibs.append(l[2:])
    #print(flibs)
    ctx.env.append_unique('STLIB', flibs)
    ctx.env.append_unique('STLIBPATH', fpaths)

def options(opt):
    opt.load('compiler_c compiler_cxx compiler_fc')
    opt.add_option('--pflotran', action='store', default=False, help='path to pflotran directory')
    opt.add_option('--crunch', action='store', default=False, help='path to crunch directory')


def setup_petsc(ctx):
   if "PETSC_DIR" not in os.environ:
      raise waflib.Errors.ConfigurationError(
         "PETSc is a required library. Please set the PETSC_DIR environment variable.")
   if "PETSC_ARCH" not in os.environ:
      raise waflib.Errors.ConfigurationError("PETSc is a required library. "
                      "Please set the PETSC_ARCH environment variable.")
   petsc_dir = os.path.abspath(os.environ["PETSC_DIR"])
   petsc_arch = os.environ["PETSC_ARCH"]
   petsc_arch_path = "{0}/{1}".format(petsc_dir, petsc_arch)
   ctx.start_msg("Using PETSc from : {0}".format(petsc_arch_path))

   # check for necessary files
   if not os.path.isdir(petsc_arch_path):
      raise waflib.Errors.ConfigurationError(
         "Could not find petsc arch directory : {0}".format(petsc_arch_path))

   petsc_variables_path = "{0}/conf/petscvariables".format(petsc_arch_path)
   if not os.path.isfile(petsc_variables_path):
      raise waflib.Errors.ConfigurationError(
         "Could not find petsc variables file : {0}".format(petsc_variables_path))

   petsc_variables = {}
   with open(petsc_variables_path, 'r') as petscfile:
      for line in petscfile:
         line = line.split("=")
         key = line[0].strip()
         value = line[1].strip()
         petsc_variables[key] = value

   if not "FC_FLAGS" in petsc_variables:
      raise waflib.Errors.ConfigurationError("PETSc : 'FC_FLAGS' missing")
   else:
      ctx.env.append_unique('FCFLAGS', petsc_variables["FC_FLAGS"].split())

   if not "PETSC_FC_INCLUDES" in petsc_variables:
      raise waflib.Errors.ConfigurationError("PETSc : 'PETSC_FC_INCLUDES' missing")
   elif petsc_variables["PETSC_FC_INCLUDES"] == "":
      raise waflib.Errors.ConfigurationError("PETSc 'PETSC_FC_INCLUDES' is empty.")
   else:
      ctx.env.append_unique('FCFLAGS', petsc_variables["PETSC_FC_INCLUDES"].split())

   if not "PETSC_CC_INCLUDES" in petsc_variables:
      raise waflib.Errors.ConfigurationError("PETSc : 'PETSC_CC_INCLUDES' missing")
   elif petsc_variables["PETSC_CC_INCLUDES"] == "":
      raise waflib.Errors.ConfigurationError("PETSc 'PETSC_CC_INCLUDES' is empty.")
   else:
      ctx.env.append_unique('CFLAGS', petsc_variables["PETSC_CC_INCLUDES"].split())
      ctx.env.append_unique('CXXFLAGS', petsc_variables["PETSC_CC_INCLUDES"].split())

   #print(ctx.env["FCFLAGS"])

   if not "PETSC_LIB_BASIC" in petsc_variables:
      raise waflib.Errors.ConfigurationError("PETSc : 'PETSC_LIB_BASIC' variable missing.")
   elif petsc_variables["PETSC_LIB_BASIC"] == "":
      raise waflib.Errors.ConfigurationError("PETSc : 'PETSC_LIB_BASIC' is empty.")
   else:
      ctx.env.append_unique("LINKFLAGS", petsc_variables["PETSC_LIB_BASIC"].split())

   if not "PETSC_EXTERNAL_LIB_BASIC" in petsc_variables:
      raise waflib.Errors.ConfigurationError("PETSc : 'PETSC_EXTERNAL_LIB_BASIC' variable missing.")
   elif petsc_variables["PETSC_EXTERNAL_LIB_BASIC"] == "":
      raise waflib.Errors.ConfigurationError("PETSc : 'PETSC_EXTERNAL_LIB_BASIC' is empty.")
   else:
      ctx.env.append_unique("LINKFLAGS", petsc_variables["PETSC_EXTERNAL_LIB_BASIC"].split())


def setup_pflotran(ctx):
   pflotran_base_dir = os.path.abspath(ctx.options.pflotran)
   pflotran_src_dir = "{0}/src/pflotran".format(pflotran_base_dir)
   pflotran_lib_short_name = "pflotranchem"
   pflotran_lib_path = "{0}/lib{1}.a".format(pflotran_src_dir,
                                             pflotran_lib_short_name)
   pflotran_definitions = "{0}/definitions.h".format(pflotran_src_dir)

   ctx.start_msg("Using pflotran from : {0}".format(pflotran_base_dir))

   # check for necessary files
   if not os.path.isfile(pflotran_lib_path):
      raise waflib.Errors.ConfigurationError(
         "Could not find pflotran chemistry library : {0}".format(pflotran_lib_path))

   if not os.path.isfile(pflotran_definitions):
      raise waflib.Errors.ConfigurationError(
         "Could not find pflotran definitions header : {0}".format(pflotran_definitions))

   ctx.env.append_unique('CXXFLAGS', ['-DHAVE_PFLOTRAN'])
   ctx.env.append_unique('STLIB', pflotran_lib_short_name)
   ctx.env.append_unique('STLIBPATH', pflotran_src_dir)
   pflotran_includes = "-I{0}".format(pflotran_src_dir)
   ctx.env.append_unique('FCFLAGS', pflotran_includes)


def setup_crunch(ctx):
        ctx.start_msg("Using crunch from : {0}".format(ctx.options.crunch))
        ctx.env.append_unique('CXXFLAGS', ['-DHAVE_CRUNCH'])
        # TODO(bja) : append -L${crunch_dir}/src/crunch -lcrunchchem to link path


def configure(ctx):
    # TODO(bja): require mpi compilers...

    # check the c compilers
    ctx.load('compiler_c')
    #ctx.check_cc()
    ctx.env['CFLAGS'] = ['-Wall', '-W', '-Wunused']

    # check the c++ compilers
    ctx.load('compiler_cxx')
    #ctx.check_cxx()
    ctx.env['CXXFLAGS'] = ['-Wall', '-W']

    # check the fortran compiler
    ctx.load('compiler_fc')
    if ctx.env.FC_NAME == 'IFORT':
        ctx.env['FCFLAGS'] = ['-warn']
    elif ctx.env.FC_NAME == 'GFORTRAN':
        ctx.env['FCFLAGS'] = ['-Wall', '-W', '-Wno-unused-parameter']

    ctx.check_fortran()
    ctx.check_fortran_verbose_flag()
    ctx.check_fortran_clib()
    ctx.check_fortran_dummy_main()
    ctx.check_fortran_mangling()

    set_fortran_stdlib(ctx)

    # both pflotran and crunch are built on petsc and (for now) we
    # require it to be present
    setup_petsc(ctx)

    if ctx.options.pflotran:
       setup_pflotran(ctx)

    if ctx.options.crunch:
       setup_crunch(ctx)


def build(ctx):
    ctx.recurse("src")

