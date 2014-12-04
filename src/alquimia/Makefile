
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

#
# Following subdirectories are built, in order.
#

SUBDIRS = \
	c \
	fortran 

# MAKECMDGOALS is the make option: make 'clobber' or 'all'
TARGET = $(MAKECMDGOALS)

# if target is undefined (i.e. MAKECMDGOALS is undefined), target = all
ifeq (,$(TARGET))
	TARGET = all
endif

#
# macro for executing TARGET in all SUBDIRS
#
ifdef SUBDIRS
$(SUBDIRS) : FORCE
	@if [ -d $@ ]; then \
		$(MAKE) --directory=$@ $(TARGET); \
	fi
endif
	@echo Build complete: $@ : $(shell date)

#
# all - Make everything in the listed sub directories
#
all : $(SUBDIRS)

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
	-rm -f .Makedepend *.o *.il

#
# FORCE - Null rule to force things to happen.
#
FORCE :

