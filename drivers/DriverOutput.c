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

#include "DriverOutput.h"
#include "alquimia/alquimia_memory.h"

typedef void (*DriverWriteMethod)(AlquimiaVectorString, AlquimiaVectorDouble*, FILE*);

struct DriverOutput 
{
  // Method for writing out file contents.
  DriverWriteMethod Write;
};

static DriverOutput* DriverOutput_New(DriverWriteMethod write_method)
{
  DriverOutput* file = malloc(sizeof(DriverOutput));
  file->Write = write_method;
  return file;
}

static void GnuplotWrite(AlquimiaVectorString var_names,
                         AlquimiaVectorDouble* var_vectors,
                         FILE* file)
{
  // Write the header.
  fprintf(file, "# ");
  for (int i = 0; i < var_names.size-1; ++i)
    fprintf(file, "\"%s\", ", var_names.data[i]);
  fprintf(file, "\"%s\"\n", var_names.data[var_names.size-1]);

  // Write the columnated data.
  for (int i = 0; i < var_vectors[0].size; ++i)
  {
    for (int j = 0; j < var_names.size; ++j)
      fprintf(file, "%14.6e", var_vectors[j].data[i]);
    fprintf(file, "\n");
  }
  fprintf(file, "\n");
}

DriverOutput* GnuplotDriverOutput_New()
{
  return DriverOutput_New(GnuplotWrite);
}

static void PythonWrite(AlquimiaVectorString var_names,
                        AlquimiaVectorDouble* var_vectors,
                        FILE* file)
{
  // Write a Python dictionary.
  fprintf(file, "data = {\n");
  for (int i = 0; i < var_names.size; ++i)
  {
    fprintf(file, "  '%s': [", var_names.data[i]);
    for (int j = 0; j < var_vectors[0].size-1; ++j)
      fprintf(file, "%14.6e, ", var_vectors[i].data[j]);
    fprintf(file, "%14.6e]", var_vectors[i].data[var_vectors[0].size-1]);
    if (i < var_names.size-1)
      fprintf(file, ",\n");
    else
      fprintf(file, "\n}\n");
  }
}

DriverOutput* PythonDriverOutput_New()
{
  return DriverOutput_New(PythonWrite);
}

void DriverOutput_WriteVectors(DriverOutput* output, 
                               const char* filename,
                               AlquimiaVectorString var_names,
                               AlquimiaVectorDouble* var_vectors)
{
  if (var_names.size == 0)
    alquimia_error("DriverOutput_WriteVectors: no variables to write!");
  for (int i = 1; i < var_names.size; ++i)
  {
    if (var_vectors[i].size != var_vectors[0].size)
    {
      alquimia_error("DriverOutput_WriteVectors: vector %d has %d data, but vector 0 has %d.", 
                     i, var_vectors[i].size, var_vectors[0].size);
    }
  }
  FILE* file = fopen(filename, "w");
  output->Write(var_names, var_vectors, file);
  fclose(file);
}

void DriverOutput_WriteMulticompVector(DriverOutput* output, 
                                       const char* filename,
                                       AlquimiaVectorString comp_names,
                                       AlquimiaVectorDouble multicomp_vector)
{
  int num_comps = comp_names.size;
  if ((multicomp_vector.size % num_comps) != 0)
  {
    alquimia_error("DriverOutput_WriteMulticompVector: multicomp_vector data has invalid size for %d components (%d).", 
                   num_comps, multicomp_vector.size);
  }
  int vec_size = multicomp_vector.size / num_comps;
  AlquimiaVectorDouble var_vectors[num_comps];
  for (int i = 0; i < num_comps; ++i)
  {
    AllocateAlquimiaVectorDouble(vec_size, &var_vectors[i]);
    for (int j = 0; j < vec_size; ++j)
      var_vectors[i].data[j] = multicomp_vector.data[num_comps*j+i];
  }
  DriverOutput_WriteVectors(output, filename, comp_names, var_vectors);
  for (int i = 0; i < num_comps; ++i)
    FreeAlquimiaVectorDouble(&var_vectors[i]);
}

