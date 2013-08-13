/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

//
// Alquimia Copyright (c) 2013, The Regents of the University of California, 
// through Lawrence Berkeley National Laboratory (subject to receipt of any 
// required approvals from the U.S. Dept. of Energy).  All rights reserved.
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

/*******************************************************************************
 **
 **  File Name: StringTokenizer.cpp
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
#include "string_tokenizer.h"

#include <sstream>
#include <string>
#include <vector>

namespace alquimia {
namespace drivers {
namespace utilities {

StringTokenizer::StringTokenizer(void) {
} /* end StringTokenizer() */

StringTokenizer::StringTokenizer(const std::string& source,
                                 const std::string& delimiters) {
  tokenize(source, delimiters);
} /* StringTokenizer(source, delimiters) */

void StringTokenizer::tokenize(const std::string& source,
                               const std::string& delimiters) {
  clear();
  std::string::size_type spos(source.find_first_not_of(delimiters, 0));
  std::string::size_type epos(source.find_first_of(delimiters, spos));

  while (std::string::npos != epos || std::string::npos != spos) {
    push_back(source.substr(spos, epos - spos));
    spos = source.find_first_not_of(delimiters, epos);
    epos = source.find_first_of(delimiters, spos);
  }
} /* end tokenize(source, delimitiers) */

void StringTokenizer::tokenize_leave_delimiters(const std::string& source,
                                               const std::string& delimiters) {
  clear();
  std::string::size_type spos(source.find_first_not_of(delimiters, 0));
  std::string::size_type epos(source.find_first_of(delimiters, spos));

  // add delimeter
  if (spos > 0) push_back(source.substr(0, spos));
  while (std::string::npos != epos || std::string::npos != spos) {
    push_back(source.substr(spos, epos - spos));
    // find position of delimiter
    spos = epos;
    epos = source.find_first_not_of(delimiters, spos);
    if (std::string::npos != epos) {
      // add delimeter
      push_back(source.substr(spos, epos - spos));
      // find position of non-delimiter
      spos = source.find_first_not_of(delimiters, epos);
      epos = source.find_first_of(delimiters, spos);
    }
  }
} /* end tokenize_with_delimiters(source, delimitiers) */

}  // namespace utilities
}  // namespace drivers
}  // namespace alquimia
