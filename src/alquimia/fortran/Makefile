# include the generic file for controling paths to dependency libraries
WORKSPACE_PATH := $(shell expr $(PWD) : '\(.*\/src\)/.*' \| $(PWD))
include $(WORKSPACE_PATH)/Makefile.inc 

HEADERS = 

SRCS = 

ifdef PFLOTRAN_DIR
	HEADERS += pflotran_alquimia.h90
	SRCS += pflotran_alquimia.F90
endif

ifdef CRUNCH_DIR
	HEADERS += crunch_alquimia.h90
	SRCS += crunch_alquimia.F90
endif

OBJS = \
	$(SRCS:%.F90=%.o)

LIBRARY = libalquimia_fortran.a

#
# Generic rules to generate object files.
#
%.o : %.F90 $(HEADERS)
	$(FC) $(FCFLAGS) -c -o $@ $<

all : $(LIBRARY)

#
# library
#
libalquimia_fortran.a : $(HEADERS) $(OBJS) $(DEPEND_LIBS)
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
	-rm -f .Makedepend *.o *.il *.mod $(LIBRARY)

#
# FORCE - Null rule to force things to happen.
#
FORCE :

