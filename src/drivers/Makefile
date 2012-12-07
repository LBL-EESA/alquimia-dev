# include the generic file for controling paths to dependency libraries
WORKSPACE_PATH := $(shell expr $(PWD) : '\(.*\/src\)/.*' \| $(PWD))
include $(WORKSPACE_PATH)/Makefile.inc 

SRCS_TEST = \
	test_main.cc

OBJS_TEST = \
	$(SRCS_TEST:%.cc=%.o)

EXECUTABLES = \
	test_main

INCLUDES = -I. \
	-I../alquimia/c++

CXXFLAGS += $(INCLUDES)

LIBS = \
	-L../alquimia/c++ -lalquimia_c++

DEPEND_LIBS += \
	../alquimia/c++/libalquimia_c++.a

#
# Generic rules to generate object files.
#
% : %.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

all : $(EXECUTABLES)

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
	rm -f .Makedepend *.o *.il $(EXECUTABLES)

#
# FORCE - Null rule to force things to happen.
#
FORCE :

