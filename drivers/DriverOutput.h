/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013-2016, The Regents of the University of California, 
** through Lawrence Berkeley National Laboratory (subject to receipt of any 
** required approvals from the U.S. Dept. of Energy).  All rights reserved.
** 
** Alquimia is available under a BSD license. See LICENSE.txt for more
** information.
**
** If you have questions about your rights to use or distribute this software, 
** please contact Berkeley Lab's Technology Transfer and Intellectual Property 
** Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
** 
** NOTICE.  This software was developed under funding from the U.S. Department 
** of Energy.  As such, the U.S. Government has been granted for itself and 
** others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
** license in the Software to reproduce, prepare derivative works, and perform 
** publicly and display publicly.  Beginning five (5) years after the date 
** permission to assert copyright is obtained from the U.S. Department of Energy, 
** and subject to any subsequent five (5) year renewals, the U.S. Government is 
** granted for itself and others acting on its behalf a paid-up, nonexclusive, 
** irrevocable, worldwide license in the Software to reproduce, prepare derivative
** works, distribute copies to the public, perform publicly and display publicly, 
** and to permit others to do so.
*/

#ifndef ALQUIMIA_DRIVER_OUTPUT_H_
#define ALQUIMIA_DRIVER_OUTPUT_H_

#include "alquimia/alquimia_containers.h"

// This type defines the interface for writing output files for Alquimia's 
// simulation drivers.
typedef struct DriverOutput DriverOutput;

// Creates a DriverOutput object that writes columnated data sensible to Gnuplot.
DriverOutput* GnuplotDriverOutput_New();

// Creates a DriverOutput object that writes a Python module that can be imported
// by an analysis script.
DriverOutput* PythonDriverOutput_New();

// Writes the given vector(s) to the file with the given name.
void DriverOutput_WriteVectors(DriverOutput* output, 
                               const char* filename,
                               AlquimiaVectorString var_names,
                               AlquimiaVectorDouble* var_vectors);

// Writes the given vector to the file, with names assigned to components 
// within the vector in a cyclic fashion.
void DriverOutput_WriteMulticompVector(DriverOutput* output, 
                                       const char* filename,
                                       AlquimiaVectorString comp_names,
                                       AlquimiaVectorDouble multicomp_vector);

#endif 
