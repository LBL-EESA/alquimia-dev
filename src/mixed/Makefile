
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

CXX = g++-mp-4.7
CXXFLAGS = -g -Wall
FC = gfortran-mp-4.7
FCFLAGS = -g -Wall

EXECUTABLES = \
	struct-derived \
	cc-f90-strings \
	struct-array-pointer

all : $(EXECUTABLES)

#
# Check compatibility of C structs and fortran derived types for POD
#
struct-derived : struct-derived.cc struct-derived.F90
	$(FC) $(FCFLAGS) -c -o $@-F90.o $@.F90
	$(CXX) $(CXXFLAGS) -c -o $@-cc.o $@.cc
	$(CXX) $(CXXFLAGS) -o $@ $@-F90.o $@-cc.o -lgfortran


#
# Pass strings (fixed length character arrays) between c++ and fortran
#
cc-f90-strings : cc-f90-strings.cc cc-f90-strings.F90
	$(FC) $(FCFLAGS) -c -o $@-F90.o $@.F90
	$(CXX) $(CXXFLAGS) -c -o $@-cc.o $@.cc
	$(CXX) $(CXXFLAGS) -o $@ $@-F90.o $@-cc.o -lgfortran

#
# Pass structs with pointers to arrays between c++ and fortran
#
struct-array-pointer : struct-array-pointer.cc struct-array-pointer.F90
	$(FC) $(FCFLAGS) -c -o $@-F90.o $@.F90
	$(CXX) $(CXXFLAGS) -c -o $@-cc.o $@.cc
	$(CXX) $(CXXFLAGS) -o $@ $@-F90.o $@-cc.o -lgfortran


#
# Passing void pointers to store fortran info
#
void-pointer : void-pointer.cc void-pointer.F90
	$(FC) $(FCFLAGS) -c -o $@-F90.o $@.F90
	$(CXX) $(CXXFLAGS) -c -o $@-cc.o $@.cc
	$(CXX) $(CXXFLAGS) -o $@ $@-F90.o $@-cc.o -lgfortran


#
# clean - Clean up the directory.
#
clean : 
	rm -f *~ *.CKP *.ln *.BAK *.bak .*.bak \
		core errs \
		,* .emacs_* \
		tags TAGS \
		make.log MakeOut

#
# clobber - Really clean up the directory.
#
clobber : clean 
	rm -f .Makedepend *.o *.il $(EXECUTABLES)

#
# FORCE - Null rule to force things to happen.
#
FORCE :

