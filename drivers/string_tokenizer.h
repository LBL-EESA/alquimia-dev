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
** 
** Authors: Benjamin Andre <bandre@lbl.gov>
*/

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
