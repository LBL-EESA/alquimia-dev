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

