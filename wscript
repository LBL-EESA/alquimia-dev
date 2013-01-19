#!/usr/bin/env python

from __future__ import print_function

import os
import string

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
    #link_flags = waflib.Tools.fc_config.parse_fortran_link(fortran_libs)
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
    #
    # check for the presence of petsc
    #
    ctx.start_msg("Checking for PETSc")
    if "PETSC_DIR" not in os.environ:
        raise waflib.Errors.ConfigurationError(
            "PETSc is a required library. "
            "Please set the PETSC_DIR environment variable.")
    if "PETSC_ARCH" not in os.environ:
        raise waflib.Errors.ConfigurationError(
            "PETSc is a required library. "
            "Please set the PETSC_ARCH environment variable.")
    petsc_dir = os.path.abspath(os.environ["PETSC_DIR"])
    petsc_arch = os.environ["PETSC_ARCH"]
    petsc_arch_path = "{0}/{1}".format(petsc_dir, petsc_arch)

    if not ctx.root.find_dir(petsc_arch_path):
        raise waflib.Errors.ConfigurationError(
            "Could not find petsc arch directory : {0}".format(petsc_arch_path))

    petsc_variables_path = "{0}/conf/petscvariables".format(petsc_arch_path)
    if not ctx.root.find_node(petsc_variables_path):
        raise waflib.Errors.ConfigurationError(
            "Could not find petsc variables file : {0}".format(petsc_variables_path))

    #
    # read the petsc variables file and store the variables as
    # key/value pairs in a dict
    #
    petsc_variables = {}
    with open(petsc_variables_path, 'r') as petscfile:
        for line in petscfile:
            line = line.split("=")
            key = line[0].strip()
            value = line[1].strip()
            petsc_variables[key] = value

    #
    # check petsc variables for fortran flags
    #
    bad_fortran = string.Template(
        "PETSc : '$FLAG' is $STATUS. Was PETSc compiled with fortran support?")
    if not "FC_FLAGS" in petsc_variables:
        raise waflib.Errors.ConfigurationError(
            bad_fortran(FLAG="FC_FLAGS", STATUS="missing"))
    else:
        ctx.env.append_unique('FCFLAGS', petsc_variables["FC_FLAGS"].split())

    if not "PETSC_FC_INCLUDES" in petsc_variables:
        raise waflib.Errors.ConfigurationError(
            bad_fortran(FLAG="PETSC_FC_INCLUDES", STATUS="missing"))
    elif petsc_variables["PETSC_FC_INCLUDES"] == "":
        raise waflib.Errors.ConfigurationError(
            bad_fortran(FLAG="PETSC_FC_INCLUDES", STATUS="empty"))
    else:
        ctx.env.append_unique('FCFLAGS', petsc_variables["PETSC_FC_INCLUDES"].split())

    #
    # check petsc variables for c flags
    #
    if not "PETSC_CC_INCLUDES" in petsc_variables:
        raise waflib.Errors.ConfigurationError("PETSc : 'PETSC_CC_INCLUDES' missing")
    elif petsc_variables["PETSC_CC_INCLUDES"] == "":
        raise waflib.Errors.ConfigurationError("PETSc 'PETSC_CC_INCLUDES' is empty.")
    else:
        ctx.env.append_unique('CFLAGS', petsc_variables["PETSC_CC_INCLUDES"].split())
        ctx.env.append_unique('CXXFLAGS', petsc_variables["PETSC_CC_INCLUDES"].split())

    #print(ctx.env["FCFLAGS"])

    #
    # check petsc variables for link libs
    #
    if not "PETSC_LIB_BASIC" in petsc_variables:
        raise waflib.Errors.ConfigurationError(
            "PETSc : 'PETSC_LIB_BASIC' variable missing.")
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
    ctx.end_msg("ok")
    ctx.msg("Using PETSc", petsc_arch_path)

def setup_pflotran(ctx):
    """
    Check for pflotran library
    """
    ctx.start_msg("Checking for PFloTran")
    pflotran_base_dir = os.path.abspath(ctx.options.pflotran)
    pflotran_src_dir = "{0}/src/pflotran".format(pflotran_base_dir)
    pflotran_lib_short_name = "pflotranchem"
    pflotran_lib_path = "{0}/lib{1}.a".format(pflotran_src_dir,
                                             pflotran_lib_short_name)
    pflotran_definitions = "{0}/definitions.h".format(pflotran_src_dir)

    # check for necessary files
    if not ctx.root.find_node(pflotran_lib_path):
        raise waflib.Errors.ConfigurationError(
            "Could not find pflotran chemistry library : {0}".format(pflotran_lib_path))

    if not ctx.root.find_node(pflotran_definitions):
        raise waflib.Errors.ConfigurationError(
            "Could not find pflotran definitions header : {0}".format(pflotran_definitions))

    ctx.env.append_unique('CFLAGS', ['-DHAVE_PFLOTRAN'])
    ctx.env.append_unique('STLIB', pflotran_lib_short_name)
    ctx.env.append_unique('STLIBPATH', pflotran_src_dir)
    pflotran_includes = "-I{0}".format(pflotran_src_dir)
    ctx.env.append_unique('FCFLAGS', pflotran_includes)
    ctx.end_msg("ok")
    ctx.msg("Using pflotran", pflotran_base_dir)
    


def setup_crunch(ctx):
    ctx.start_msg("Using crunch from : {0}".format(ctx.options.crunch))
    ctx.env.append_unique('CFLAGS', ['-DHAVE_CRUNCH'])
    # TODO(bja) : append -L${crunch_dir}/src/crunch -lcrunchchem to link path


def configure(ctx):
    ctx.define("ALQUIMIA_MAJOR_VERSION", '0')
    ctx.define("ALQUIMIA_MINOR_VERSION", '1')
    ctx.define("ALQUIMIA_PATCH_VERSION", '0')

    # TODO(bja): require mpi compilers...

    # TODO(bja): use debug by default, turn off via cmdline flag

    #
    # check the c compilers
    #
    ctx.load('compiler_c')
    #ctx.check_cc()
    ctx.check_cc(header_name='stdlib.h')
    ctx.check_cc(function_name='calloc', header_name='stdlib.h')
    ctx.check_cc(function_name='free', header_name='stdlib.h')
    ctx.check_cc(header_name='stdio.h')
    ctx.check_cc(function_name='fprintf', header_name='stdio.h')
    ctx.check_cc(header_name='string.h')
    ctx.check_cc(function_name='strncpy', header_name='string.h')
    ctx.check_cc(function_name='strlen', header_name='string.h')
    ctx.check_cc(header_name='ctype.h')
    ctx.check_cc(header_name='mpi.h')
    ctx.check_cc(lib='mpi')
    ctx.check_cc(function_name='MPI_Init', header_name='mpi.h')
    ctx.env['CFLAGS'] = ['-Wall', '-W', '-Wunused', '-g']

    #
    # check the c++ compilers
    #
    ctx.load('compiler_cxx')
    #ctx.check_cxx()
    ctx.check_cxx(header_name='cstdlib')
    ctx.check_cxx(header_name='cctype')
    ctx.check_cxx(header_name='cstring')
    ctx.check_cxx(header_name='unistd.h', manditory=False)
    ctx.check_cxx(header_name='string')
    ctx.check_cxx(header_name='sstream')
    ctx.check_cxx(header_name='fstream')
    ctx.check_cxx(header_name='iostream')
    ctx.check_cxx(header_name='iomanip')
    ctx.check_cxx(header_name='vector')
    ctx.check_cxx(header_name='map')
    ctx.check_cxx(header_name='stdexcept')
    ctx.check_cxx(header_name='mpi.h')
    ctx.check_cxx(lib='mpi')
    ctx.check_cxx(function_name='MPI_Init', header_name='mpi.h')
    #ctx.check_cxx(header_name='')
    ctx.env['CXXFLAGS'] = ['-Wall', '-W', '-g']

    #
    # check the fortran compiler
    #
    ctx.load('compiler_fc')
    if ctx.env.FC_NAME == 'IFORT':
        ctx.env['FCFLAGS'] = ['-warn']
    elif ctx.env.FC_NAME == 'GFORTRAN':
        ctx.env['FCFLAGS'] = ['-Wall', '-W', '-Wno-unused-parameter', '-g']
        #ctx.env['FCFLAGS'] = ['-std=f2003', '-Wall', '-W', '-Wno-unused-parameter', '-g']

    ctx.check_fortran()
    ctx.check_fortran_verbose_flag()
    ctx.check_fortran_dummy_main()
    ctx.check_fortran_mangling()
    set_fortran_stdlib(ctx)

    #
    # engine libraries
    #

    # both pflotran and crunch are built on petsc and (for now) we
    # require it to be present
    setup_petsc(ctx)

    if ctx.options.pflotran:
       setup_pflotran(ctx)

    if ctx.options.crunch:
       setup_crunch(ctx)

    #
    # write config header and linking file
    #

    ctx.write_config_header('alquimia_config.h')


def build(ctx):
    ctx.recurse("src")

    # create the variables file
    bld_dir = ctx.path.get_bld()
    aq_vars = bld_dir.make_node("alquimia_variables")
    ctx(rule='touch ${TGT}', target=aq_vars)
    # write variable info
    aq_vars.write("# include and link flags for compiling against alquimia and the available engines.\n", flags='w')
    aq_vars.write("CFLAGS = -I{0}/include\n".format(ctx.env.PREFIX), flags='a')
    aq_vars.write("CXXFLAGS = -I{0}/include\n".format(ctx.env.PREFIX), flags='a')
    # TODO(bja) : list of link libs
    link_flags = ""
    aq_vars.write("LIBS = {0}\n".format(link_flags), flags='a')
    # install the variables file
    ctx.install_files("${PREFIX}", "alquimia_variables")
