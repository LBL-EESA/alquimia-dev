/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <limits>

struct data {
  int xlength;
  float* x;
  int ylength;
  double* y;
 };


extern "C" {
  void process_data_ (data* p) ;
}

void print_data(const data& d) {
  std::cout << "x_length : " << d.xlength << std::endl;
  for (int i = 0; i < d.xlength; ++i) {
    std::cout << d.x[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "y_length : " << d.ylength << std::endl;
  for (int i = 0; i < d.ylength; ++i) {
    std::cout << d.y[i] << " ";
  }
  std::cout << std::endl;
}

int main(int argc, char* argv[])
{
  int length = 10;
  if (argc > 1) {
    std::stringstream temp(argv[1]);
    temp >> length;
  }
  std::cout << "Using length : " << length << std::endl;
  data d;
  d.xlength = length;
  d.x = new float[d.xlength];
  for (int i = 0; i < d.xlength; ++i) {
    d.x[i] = i;
  }
  d.ylength = 0;
  d.y = NULL;
  print_data(d);
  process_data_(&d);
  print_data(d);
  return EXIT_SUCCESS;
}
