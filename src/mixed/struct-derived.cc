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
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <limits>

struct point {
  long w;
  float x;
  int y;
  double z;
 };

void print_point(const point& p) {
  std::cout << "(w, x, y, z) = ("
            << p.w << ", " << p.x << ", " << p.y << ", " << p.z << ")\n";
}

extern "C" {
  void negate_ (point* p) ;
}

int main(int argc, char* argv[])
{
  point p;
  p.w = std::numeric_limits<long>::max();
  p.x = std::numeric_limits<float>::epsilon();
  p.y = std::numeric_limits<int>::min();
  p.z = std::numeric_limits<double>::max();
  print_point(p);
  negate_(&p);
  print_point(p);
  return EXIT_SUCCESS;
}

