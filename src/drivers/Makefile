
# Alquimia Copyright (c) 2013, The Regents of the University of California, 
# through Lawrence Berkeley National Laboratory (subject to receipt of any 
# required approvals from the U.S. Dept. of Energy).  All rights reserved.
# 
# Alquimia is available under a BSD license. See LICENSE.txt for more
# information.
#
# If you have questions about your rights to use or distribute this software, 
# please contact Berkeley Lab's Technology Transfer and Intellectual Property 
# Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
# 
# NOTICE.  This software was developed under funding from the U.S. Department 
# of Energy.  As such, the U.S. Government has been granted for itself and 
# others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
# license in the Software to reproduce, prepare derivative works, and perform 
# publicly and display publicly.  Beginning five (5) years after the date 
# permission to assert copyright is obtained from the U.S. Department of Energy, 
# and subject to any subsequent five (5) year renewals, the U.S. Government is 
# granted for itself and others acting on its behalf a paid-up, nonexclusive, 
# irrevocable, worldwide license in the Software to reproduce, prepare derivative
# works, distribute copies to the public, perform publicly and display publicly, 
# and to permit others to do so.
# 
# Authors: Benjamin Andre <bandre@lbl.gov>

# include the generic file for controlling paths to dependency libraries
WORKSPACE_PATH := $(shell expr $(PWD) : '\(.*\/alquimia\/src\/drivers\)/.*' \| $(PWD))
include $(WORKSPACE_PATH)/Makefile.inc 

# TODO(bja): cleanup dependencies using something like makedepend

# demo drives = always in debug mode!
CXXFLAGS += -g
# TODO(bja): need to make this generic!
#LDFLAGS +=  -L/opt/local/lib/gcc46 -lgfortran
# NOTE(bja): seems to work without specifying the search path...
LDFLAGS +=  -lgfortran

#
# generic to include alquimia library
#
INCLUDES = -I. \
	-I../alquimia/c

CXXFLAGS += $(INCLUDES)

LIBS = \
	-L../alquimia/c -lalquimia_c \
	-L../alquimia/c -lalquimia_cutils \
	-L../alquimia/fortran -lalquimia_fortran

DEPEND_LIBS += \
	../alquimia/fortran/libalquimia_fortran.a \
	../alquimia/c/libalquimia_c.a \
	../alquimia/c/libalquimia_cutils.a

#
# libdemoutils.a
#
SRCS_DEMOUTIL = \
	demo_utils.cc \
	demo_containers.cc \
	demo_output.cc \
	cfg_reader.cc \
	string_tokenizer.cc

OBJS_DEMOUTIL = \
	$(SRCS_DEMOUTIL:%.cc=%.o)

HEADERS_DEMOUTIL = \
	demo_utils.h \
	demo_containers.h \
	demo_output.h \
	cfg_reader.h \
	string_tokenizer.h

LIBRARIES_DEMOUTIL = libdemoutils.a

#
# C++ batch_chem demo driver
#
SRCS_BATCH = \
	batch_chem.cc

OBJS_BATCH = \
	$(SRCS_BATCH:%.cc=%.o)

HEADERS_BATCH = \
	batch_chem.h

EXECUTABLES_BATCH = \
	batch_chem

LIBS_BATCH = $(LIBS) \
	-L. -ldemoutils

DEPEND_LIBS_BATCH = $(DEPEND_LIBS) \
	$(LIBRARIES_DEMOUTIL)

#
# test_alquimia_interface_factory
#
SRCS_TEST = \
	test_alquimia_interface_factory.cc

OBJS_TEST = \
	$(SRCS_TEST:%.cc=%.o)

EXECUTABLES_TEST = \
	test_alquimia_interface_factory

#
# Generic rules to generate object files.
#
%.o : %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

all : batch_chem

#
# demo utils library
#
libdemoutils.a : $(HEADERS_DEMOUTIL) $(OBJS_DEMOUTIL)
	$(AR) $(ARFLAGS) $@ $(OBJS_DEMOUTIL)

#
# batch_chem_cc
#
batch_chem : $(HEADERS_BATCH) $(OBJS_BATCH) $(DEPEND_LIBS_BATCH)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_BATCH) $(LIBS_BATCH) $(LDFLAGS)


#
# clean - Clean up the directory.
#
clean : $(SUBDIRS)
	rm -f *~ *.CKP *.ln *.BAK *.bak .*.bak \
		core errs \
		,* .emacs_* \
		tags TAGS \
		make.log MakeOut

#
# clobber - Really clean up the directory.
#
clobber : clean $(SUBDIRS)
	rm -f .Makedepend *.o *.il \
		$(EXECUTABLES_TEST) $(EXECUTABLES_BATCH) \
		$(LIBRARIES_DEMOUTIL)

#
# FORCE - Null rule to force things to happen.
#
FORCE :

