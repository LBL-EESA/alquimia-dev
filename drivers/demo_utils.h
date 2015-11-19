/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013-2015, The Regents of the University of California, 
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
** 
** Authors: Benjamin Andre <bandre@lbl.gov>
*/

#ifndef ALQUIMIA_DRIVERS_CCDEMOUTILS_H_
#define ALQUIMIA_DRIVERS_CCDEMOUTILS_H_

#include <vector>
#include <map>
#include <string>
#include <iomanip>
#include <iostream>

//
// Common stand alone utility functions
//

namespace alquimia {
namespace drivers {
namespace utilities {

/*******************************************************************************
 **
 **  Custom comparison operators
 **
 ******************************************************************************/
bool CaseInsensitiveStringCompare(const std::string& string1, 
                                  const std::string& string2);
bool CompareFabs(const double& a, const double& b);

/*******************************************************************************
 **
 **  String conversion utilities
 **
 ******************************************************************************/
void LowerCaseString(const std::string& in, std::string* out);
void RemoveLeadingAndTrailingWhitespace(std::string* line);

/*******************************************************************************
 **
 **  Print Utilities
 **
 ******************************************************************************/
template <typename T>
void PrintVector(const std::string& name, 
                 const std::vector<T>& data,
                 const int precision = -1,
                 const bool comma_seperated = false) {
  if (precision > 0) {
    std::cout << std::setprecision(precision);
  }
  std::cout << name << " : { ";
  for (typename std::vector<T>::const_iterator i = data.begin();
       i != data.end(); ++i) {
    std::cout << *i;
    if (i != --data.end()) {
      if (comma_seperated) {
        std::cout << ", ";
      } else {
        std::cout << "  ";
      }
    }
  }
  std::cout << " }\n";
}  // end PrintVector

template <typename T>
void PrintMap(const std::string& name, 
              const std::map<std::string, T>& data,
              const int precision = -1,
              const bool comma_seperated = false) {
  if (precision > 0) {
    std::cout << std::setprecision(precision);
  }
  std::cout << name << " : { ";
  for (typename std::map<std::string, T>::const_iterator i = data.begin();
       i != data.end(); ++i) {
    std::cout << "'" << i->first << "' : " << i->second;
    if (i != --data.end()) {
      if (comma_seperated) {
        std::cout << ", ";
      } else {
        std::cout << "  ";
      }
    }
  }
  std::cout << " }\n";
}  // end PrintVector


}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
#endif  // ALQUIMIA_DRIVERS_CCDEMOUTILS_H_
