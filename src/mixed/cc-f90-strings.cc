/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

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
