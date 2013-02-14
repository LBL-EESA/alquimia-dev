/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_DRIVERS_CCDEMOUTILS_H_
#define ALQUIMIA_DRIVERS_CCDEMOUTILS_H_

#include <vector>
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


}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
#endif  // ALQUIMIA_DRIVERS_CCDEMOUTILS_H_
