# include the generic file for controling paths to dependency libraries
WORKSPACE_PATH := $(shell expr $(PWD) : '\(.*\/src\)/.*' \| $(PWD))
include $(WORKSPACE_PATH)/Makefile.inc 

INCLUDES = -I../alquimia/c -I.

CFLAGS += $(INCLUDES)

LIBS = \
	-L../alquimia/c -lalquimia_c \
	-L../alquimia/c -lalquimia_cutils \
	-L../alquimia/fortran -lalquimia_fortran

DEPEND_LIBS += \
	../alquimia/fortran/libalquimia_fortran.a \
	../alquimia/c/libalquimia_c.a \
	../alquimia/c/libalquimia_cutils.a

#
# test for alquimia c utilities
#

UTIL_TEST = test_alquimia_c_utils

UTIL_TEST_SRCS = test_alquimia_c_utils.c

UTIL_TEST_OBJS = \
	$(UTIL_TEST_SRCS:%.c=%.o)

#
# Generic rules to generate object files.
#
%.o : %.c $(HEADERS) $(UTIL_HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<


all : $(UTIL_TEST)

test_alquimia_c_utils : $(UTIL_TEST_OBJS) $(DEPEND_LIBS)
	$(CC) $(CFLAGS) -o $@ $(UTIL_TEST_OBJS) $(UTIL_LIBRARY) $(LIBS) $(LDFLAGS)

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

