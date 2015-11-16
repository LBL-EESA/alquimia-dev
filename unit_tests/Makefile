
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

