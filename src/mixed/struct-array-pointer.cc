/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

//
// Alquimia Copyright (c) 2013, The Regents of the University of California, 
// through Lawrence Berkeley National Laboratory (subject to receipt of any 
// required approvals from the U.S. Dept. of Energy).  All rights reserved.
// 
// Alquimia is available under a BSD license. See LICENSE.txt for more
// information.
//
// If you have questions about your rights to use or distribute this software, 
// please contact Berkeley Lab's Technology Transfer and Intellectual Property 
// Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
// 
// NOTICE.  This software was developed under funding from the U.S. Department 
// of Energy.  As such, the U.S. Government has been granted for itself and 
// others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
// license in the Software to reproduce, prepare derivative works, and perform 
// publicly and display publicly.  Beginning five (5) years after the date 
// permission to assert copyright is obtained from the U.S. Department of Energy, 
// and subject to any subsequent five (5) year renewals, the U.S. Government is 
// granted for itself and others acting on its behalf a paid-up, nonexclusive, 
// irrevocable, worldwide license in the Software to reproduce, prepare derivative
// works, distribute copies to the public, perform publicly and display publicly, 
// and to permit others to do so.
// 
// Authors: Benjamin Andre <bandre@lbl.gov>
//


#include <cstdlib>
#include <cstring>

#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <limits>

#define MAX_STRING_LENGTH 256

struct container {
  int num;
  struct data* d;
};

struct data {
  char* name;
  char* outstr;
  int xlength;
  int* x;
  int ylength;
  double* y;
 };


extern "C" {
  void process_container_ (container* p) ;
}

void initialize_container(container* c, int num, int length) {
  c->num = num;
  c->d = new data[c->num];

  for (int i = 0; i < c->num; ++i) {
    data* d = &(c->d[i]);
    std::stringstream name;
    name << "data_" << i;
    d->name = new char[MAX_STRING_LENGTH];
    strncpy(d->name, name.str().c_str(), MAX_STRING_LENGTH);

    d->outstr = new char[MAX_STRING_LENGTH];

    d->xlength = length - num + i;
    d->x = new int[d->xlength];
    for (int j = 0; j < d->xlength; ++j) {
      d->x[j] = 2 * j - i;
    }

    d->ylength = 0;
    d->y = NULL;
  }
}

void free_container(container* c) {
  for (int i = 0; i < c->num; ++i) {
    delete c->d[i].name;
    delete c->d[i].outstr;
    delete c->d[i].x;
    delete c->d[i].y;
  }
}

void print_container(const container& c) {
  std::cout << " container->num : " << c.num << std::endl;
  for (int i = 0; i < c.num; ++i) {
    struct data* d = &(c.d[i]);
    std::cout << "  d->name : " << d->name << std::endl;
    if (d->outstr == NULL) {
      std::cout << "  d->outstr : NULL" << std::endl;
    } else {
      std::cout << "  d->outstr : '" << d->outstr << "'" << std::endl;
    }
    std::cout << "    d->x_length : " << d->xlength << "\n     ";
    for (int j = 0; j < d->xlength; ++j) {
      std::cout << d->x[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "    d->y_length : " << d->ylength << "\n     ";
    for (int j = 0; j < d->ylength; ++j) {
      std::cout << d->y[j] << " ";
    }
    std::cout << std::endl;
  }
}

int main(int argc, char* argv[])
{
  int num = 2;  // number of data structs
  int length = 10;  // length of arrays
  if (argc == 3) {
    std::stringstream temp1(argv[1]);
    temp1 >> num;
    std::stringstream temp2(argv[2]);
    temp2 >> length;
  }
  std::cout << "Using num : " << num << std::endl;
  std::cout << "Using length : " << length << std::endl;

  container c;
  initialize_container(&c, num, length);

  std::cout << "\nbefore fortran :\n";
  print_container(c);

  process_container_(&c);

  std::cout << "\nafter fortran :\n";
  print_container(c);

  free_container(&c);

  return EXIT_SUCCESS;
}
