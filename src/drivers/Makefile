# include the generic file for controling paths to dependency libraries
WORKSPACE_PATH := $(shell expr $(PWD) : '\(.*\/src\)/.*' \| $(PWD))
include $(WORKSPACE_PATH)/Makefile.inc 

# TODO(bja): cleanup dependencies using something like makedepend

# demo drives = always in debug mode!
CXXFLAGS += -g
# TODO(bja): need to make this generic!
LDFLAGS +=  -L/opt/local/lib/gcc46 -lgfortran

#
# generic to include alquimia library
#
INCLUDES = -I. \
	-I../alquimia/c \
	-I../alquimia/c++

CXXFLAGS += $(INCLUDES)

LIBS = \
	-L../alquimia/c++ -lalquimia_c++ \
	-L../alquimia/fortran -lalquimia_fortran \
	-L../alquimia/c -lalquimia_c

DEPEND_LIBS += \
	../alquimia/c++/libalquimia_c++.a \
	../alquimia/fortran/libalquimia_fortran.a \
	../alquimia/c/libalquimia_c.a 

#
# libdemoutils.a
#
SRCS_DEMOUTIL = \
	demo_utils.cc \
	demo_containers.cc \
	cfg_reader.cc \
	string_tokenizer.cc

OBJS_DEMOUTIL = \
	$(SRCS_DEMOUTIL:%.cc=%.o)

HEADERS_DEMOUTIL = \
	demo_utils.h \
	demo_containers.h \
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

all : test_alquimia_interface_factory batch_chem

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
# test program
#
test_alquimia_interface_factory : $(OBJS_TEST) $(DEPEND_LIBS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_TEST) $(LIBS) $(LDFLAGS)

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

