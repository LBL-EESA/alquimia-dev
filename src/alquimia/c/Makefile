# include the generic file for controling paths to dependency libraries
WORKSPACE_PATH := $(shell expr $(PWD) : '\(.*\/src\)/.*' \| $(PWD))
include $(WORKSPACE_PATH)/Makefile.inc 

INCLUDES = -I.

CFLAGS += $(INCLUDES)

#
# alquimia c library
#
HEADERS = \
	alquimia_constants.h \
	alquimia_containers.h \
	pflotran_alquimia_interface.h

SRCS = \
	alquimia_constants.c 


OBJS = \
	$(SRCS:%.c=%.o)

LIBRARY = libalquimia_c.a

#
# alquimia c utilities library
#
UTIL_HEADERS = \
	alquimia_interface.h \
	alquimia_memory.h \
	alquimia_util.h 

UTIL_SRCS = \
	alquimia_interface.c \
	alquimia_memory.c \
	alquimia_util.c

UTIL_OBJS = \
	$(UTIL_SRCS:%.c=%.o)

UTIL_LIBRARY = libalquimia_cutils.a

#
# Generic rules to generate object files.
#
%.o : %.c $(HEADERS) $(UTIL_HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<


all : $(LIBRARY) $(UTIL_LIBRARY) $(UTIL_TEST)

libalquimia_c.a : $(HEADERS) $(OBJS) $(DEPEND_LIBS)
	$(AR) $(ARFLAGS) $@ $(OBJS)


libalquimia_cutils.a : $(HEADERS) $(UTIL_HEADERS) $(UTIL_OBJS) $(DEPEND_LIBS)
	$(AR) $(ARFLAGS) $@ $(UTIL_OBJS)

#
# clean - Clean up the directory.
#
clean :
	-rm -f *~ *.CKP *.ln *.BAK *.bak .*.bak \
		core errs \
		,* .emacs_* \
		tags TAGS \
		make.log MakeOut

#
# clobber - Really clean up the directory.
#
clobber : clean
	-rm -f .Makedepend *.o *.il *.mod $(LIBRARY) $(UTIL_LIBRARY)

#
# FORCE - Null rule to force things to happen.
#
FORCE :

