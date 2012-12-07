# include the generic file for controling paths to dependency libraries
WORKSPACE_PATH := $(shell expr $(PWD) : '\(.*\/src\)/.*' \| $(PWD))
include $(WORKSPACE_PATH)/Makefile.inc 

HEADERS = \
	alquimia_containers.h \
	alquimia_strings.h \
	alquimia_interface_factory.h \
	alquimia_interface.h

SRCS = \
	alquimia_interface_factory.cc 

ifdef PFLOTRAN_DIR
	HEADERS += pflotran_alquimia_interface.h
	SRCS += pflotran_alquimia_interface.cc
endif

ifdef CRUNCH_DIR
	HEADERS += crunch_alquimia_interface.h
	SRCS += crunch_alquimia_interface.cc 
endif

OBJS = \
	$(SRCS:%.cc=%.o)

LIBRARY = libalquimia_c++.a

INCLUDES = -I.

CXXFLAGS += $(INCLUDES)

#
# Generic rules to generate object files.
#
% : %.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $<

all : $(LIBRARY)

#
# library
#
libalquimia_c++.a : $(HEADERS) $(OBJS) $(DEPEND_LIBS)
	$(AR) $(ARFLAGS) $@ $(OBJS)


#
# clean - Clean up the directory.
#
clean : $(SUBDIRS)
	-rm -f *~ *.CKP *.ln *.BAK *.bak .*.bak \
		core errs \
		,* .emacs_* \
		tags TAGS \
		make.log MakeOut

#
# clobber - Really clean up the directory.
#
clobber : clean $(SUBDIRS)
	-rm -f .Makedepend *.o *.il $(LIBRARY)

#
# FORCE - Null rule to force things to happen.
#
FORCE :

