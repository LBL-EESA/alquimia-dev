/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "cc_demo_utils.h"

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
    if (std::tolower(*char1) != std::tolower(*char2)) {
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
    *c = std::tolower(*c);
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
