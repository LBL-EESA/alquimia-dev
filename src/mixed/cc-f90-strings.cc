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
#include <string>
#include <cstring>

#include "cc-f90-strings.h"

extern "C" {
  void reverse_string_ (char in_str[STRING_MAX_LENGTH],
                        int& in_length,
                        char out_str[STRING_MAX_LENGTH]) ;
}

void clear_string(char str[]) {
  for (int i = 0; i < STRING_MAX_LENGTH; ++i) {
    str[i] = '\0';
  }
}

int main(int argc, char* argv[])
{
  std::string filename("data-file.in");
  char in_str[STRING_MAX_LENGTH];
  char out_str[STRING_MAX_LENGTH];

  // NOTE: it appears to be necessary to nullify the entire string
  // prior to passing them to fortran....
  //clear_string(in_str);
  //clear_string(out_str);

  int length = filename.copy(in_str, filename.length(), 0);
  in_str[length] = '\0';
  std::cout << "length : " << length << std::endl;
  std::cout << "before:\n  '" << in_str << "'\n  '" << out_str << "'" << std::endl;

  reverse_string_(in_str, length, out_str);

  std::cout << "after:\n  '" << in_str << "'\n  '" << out_str << "'" << std::endl;
  
  if (false) {
    for (int i = 0; i < STRING_MAX_LENGTH; ++i) {
      std::cout << "'" << out_str[i];
    }
  }
  return EXIT_SUCCESS;
}
