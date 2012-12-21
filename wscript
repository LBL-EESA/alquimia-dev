#!/usr/bin/env python

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
    print(flibs)
    ctx.env.append_unique('STLIB', flibs)
    ctx.env.append_unique('STLIBPATH', fpaths)

def options(opt):
    opt.load('compiler_c compiler_cxx compiler_fc')
    opt.add_option('--pflotran', action='store', default=False, help='path to pflotran directory')
    opt.add_option('--crunch', action='store', default=False, help='path to crunch directory')

def configure(ctx):

    # check the c compilers
    ctx.load('compiler_c')
    #ctx.check_cc()
    ctx.env['CFLAGS'] = ['-Wall', '-W']

    # check the c++ compilers
    ctx.load('compiler_cxx')
    #ctx.check_cxx()
    ctx.env['CXXFLAGS'] = ['-Wall', '-W']

    # check the fortran compiler
    ctx.load('compiler_fc')
    if ctx.env.FC_NAME == 'IFORT':
        ctx.env['FCFLAGS'] = ['-warn']
    elif ctx.env.FC_NAME == 'GFORTRAN':
        ctx.env['FCFLAGS'] = ['-Wall', '-W']

    ctx.check_fortran()
    ctx.check_fortran_verbose_flag()
    ctx.check_fortran_clib()
    ctx.check_fortran_dummy_main()
    ctx.check_fortran_mangling()

    set_fortran_stdlib(ctx)

    if ctx.options.pflotran:
        print "--- using pflotran dir = {0}".format(ctx.options.pflotran)
        ctx.env.append_unique('CXXFLAGS', ['-DHAVE_PFLOTRAN'])
        # TODO(bja) : append -L${pflotran_dir}/src/pflotran -lpflotranchem to link path

    if ctx.options.crunch:
        print "--- using crunch dir = {0}".format(ctx.options.crunch)
        ctx.env.append_unique('CXXFLAGS', ['-DHAVE_CRUNCH'])
        # TODO(bja) : append -L${crunch_dir}/src/crunch -lcrunchchem to link path

def build(ctx):
    ctx.recurse("src")

