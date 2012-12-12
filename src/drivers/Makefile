# include the generic file for controling paths to dependency libraries
WORKSPACE_PATH := $(shell expr $(PWD) : '\(.*\/src\)/.*' \| $(PWD))
include $(WORKSPACE_PATH)/Makefile.inc 

# TODO(bja): cleanup dependencies using something like makedepend

# demo drives = always in debug mode!
CXXFLAGS += -g

#
# generic to include alquimia library
#
INCLUDES = -I. \
	-I../alquimia/c++

CXXFLAGS += $(INCLUDES)

LIBS = \
	-L../alquimia/c++ -lalquimia_c++

DEPEND_LIBS += \
	../alquimia/c++/libalquimia_c++.a

#
# libdemoutils_c++.a
#
SRCS_CCDEMOUTIL = \
	cc_demo_utils.cc \
	cfg_reader.cc \
	string_tokenizer.cc

OBJS_CCDEMOUTIL = \
	$(SRCS_CCDEMOUTIL:%.cc=%.o)

HEADERS_CCDEMOUTIL = \
	cc_demo_utils.h \
	cfg_reader.h \
	string_tokenizer.h

LIBRARIES_CCDEMOUTIL = libdemoutils_c++.a

#
# C++ batch_chem demo driver
#
SRCS_CCBATCH = \
	cc_batch_chem.cc

OBJS_CCBATCH = \
	$(SRCS_CCBATCH:%.cc=%.o)

HEADERS_CCBATCH = \
	cc_batch_chem.h

EXECUTABLES_CCBATCH = \
	cc_batch_chem

LIBS_CCBATCH = $(LIBS) \
	-L. -ldemoutils_c++

DEPEND_LIBS_CCBATCH = $(DEPEND_LIBS) \
	$(LIBRARIES_CCDEMOUTIL)

#
# test_main
#
SRCS_TEST = \
	test_main.cc

OBJS_TEST = \
	$(SRCS_TEST:%.cc=%.o)

EXECUTABLES_TEST = \
	test_main

#
# Generic rules to generate object files.
#
% : %.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

all : test_main cc_batch_chem

#
# demo utils library
#
libdemoutils_c++.a : $(HEADERS_CCDEMOUTIL) $(OBJS_CCDEMOUTIL)
	$(AR) $(ARFLAGS) $@ $(OBJS_CCDEMOUTIL)

#
# batch_chem_cc
#
cc_batch_chem : $(HEADERS_CCBATCH) $(OBJS_CCBATCH) $(DEPEND_LIBS_CCBATCH)
	echo $(HEADERS_CCBATCH) $(OBJS_CCBATCH) $(DEPEND_LIBS_CCBATCH)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_CCBATCH) $(LIBS_CCBATCH) $(LDFLAGS)

#
# test program
#
test_main : $(OBJS_TEST) $(DEPEND_LIBS)
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
		$(EXECUTABLES_TEST) $(EXECUTABLES_CCBATCH) \
		$(LIBRARIES_CCDEMOUTIL)

#
# FORCE - Null rule to force things to happen.
#
FORCE :

