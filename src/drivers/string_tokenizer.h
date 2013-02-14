/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef ALQUIMIA_DRIVERS_STRING_TOKENIZER_HH_
#define ALQUIMIA_DRIVERS_STRING_TOKENIZER_HH_

/*******************************************************************************
 **
 **  File Name: StringTokenizer.h
 **
 **  Source Control: $Id$
 **
 **  Description: This class is a simple string tokenizer.
 **
 **  Notes:
 **    - Found July 2003, unencumbered, on the web at:
 **
 **      http:// www.thecodezone.com/diary/archives/000057.html
 **
 **    - no longer available from origional source, but archived at:
 **
 **      http:// web.archive.org/web/20030810163805/http:// www.thecodezone.com/diary/archives/000057.html
 **
 *******************************************************************************/

#include <string>
#include <vector>

namespace alquimia {
namespace drivers {
namespace utilities {

class StringTokenizer : public std::vector<std::string> {
 public:

  StringTokenizer(void);
  StringTokenizer(const std::string& source,
                  const std::string& delimiters = " \t\n");
  void tokenize(const std::string& source,
                const std::string& delimiters = " \t\n");
  // the following tokenizes, but places delimiters in list too - geh
  void tokenize_leave_delimiters(const std::string& source,
                                const std::string& delimiters = " \t\n");
};

}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
#endif     /* ALQUIMIA_DRIVERS_STRING_TOKENIZER_HH_ */
