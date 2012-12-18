# include the generic file for controling paths to dependency libraries
WORKSPACE_PATH := $(shell expr $(PWD) : '\(.*\/src\)/.*' \| $(PWD))
include $(WORKSPACE_PATH)/Makefile.inc 

HEADERS = \
	alquimia_containers.h \
	alquimia_memory.h \
	alquimia_util.h

SRCS = \
	alquimia_memory.c \
	alquimia_util.c

OBJS = \
	$(SRCS:%.c=%.o)

LIBRARY = libalquimia_c.a

INCLUDES = -I../ -I.

CXXFLAGS += $(INCLUDES)

#
# Generic rules to generate object files.
#
%.o : %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<


#
# library
#
all : $(LIBRARY)

libalquimia_c.a : $(HEADERS) $(OBJS) $(DEPEND_LIBS)
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

