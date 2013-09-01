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

#include "demo_utils.h"

#include <cmath>

#include <string>
#include <sstream>
#include <vector>

namespace alquimia {
namespace drivers {
namespace utilities {

/*******************************************************************************
 **
 **  Custom comparison operators
 **
 ******************************************************************************/
bool CaseInsensitiveStringCompare(const std::string& string1, 
                                  const std::string& string2) {
  // really CaseInsensitiveStringsEqual...

  // if not the same length, not the same. don't bother checking
  // individual characters
  if (string1.size() != string2.size()) {
    return false;
  }
  // loop through each string, check each character individually.
  std::string::const_iterator char1, char2;
  for (char1 = string1.begin(), char2 = string2.begin(); 
       char1 != string1.end(); ++char1, ++char2) {
    if (tolower(*char1) != tolower(*char2)) {
      return false;
    }
  }
  return true;
}  // end CaseInsensitiveStringCompare()

      // std::stringstream output;
      // output << "CICS : strings are not the same length\n"
      //        << "  string one : " << string1
      //        << "\n  string two : " << string2 << std::endl;
      // chem_out.Write(kVerbose, output.str());

      // std::stringstream output;
      // output << "CICS : strings differ at position "
      //        << std::distance(string1.begin(), char1)
      //        << "\n  string one value: " << *char1
      //        << "\n  string two value: " << *char2 << std::endl;
      // chem_out.Write(kVerbose, output.str());


bool CompareFabs(const double& a, const double& b) {
  /* for use with stl algorithms such as max element. Must "return
   * true if the first argument is to be considered less than the
   * second argument, and false otherwise" */
  return std::fabs(a) < std::fabs(b);
}  // end CompareFabs()


/*******************************************************************************
 **
 **  Conversion routines
 **
 ******************************************************************************/

/*
** convert string to lower case
*/
void LowerCaseString(const std::string& in, std::string* out) {
  *out = in;
  for (std::string::iterator c = out->begin(); c != out->end(); ++c) {
    *c = tolower(*c);
  }
}

void RemoveLeadingAndTrailingWhitespace(std::string* line) {
  std::string whitespace(" \t\f\v\n\r");
  size_t start = line->find_first_not_of(whitespace);
  if (start != std::string::npos) {
    line->erase(0, start);
  } else if (start == std::string::npos) {
    // entire line is blank
    line->erase(0);
  }
  size_t end = line->find_last_not_of(whitespace);
  if (end != std::string::npos) {
    ++end;  // find returned the last non-whitespace character....
    line->erase(end);
  }
}  // end RemoveLeadingAndTrailingWhitespace()

}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
