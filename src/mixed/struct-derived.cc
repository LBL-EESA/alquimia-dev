/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <limits>

struct point {
  long w;
  float x;
  int y;
  double z;
 };

void print_point(const point& p) {
  std::cout << "(w, x, y, z) = ("
            << p.w << ", " << p.x << ", " << p.y << ", " << p.z << ")\n";
}

extern "C" {
  void negate_ (point* p) ;
}

int main(int argc, char* argv[])
{
  point p;
  p.w = std::numeric_limits<long>::max();
  p.x = std::numeric_limits<float>::epsilon();
  p.y = std::numeric_limits<int>::min();
  p.z = std::numeric_limits<double>::max();
  print_point(p);
  negate_(&p);
  print_point(p);
  return EXIT_SUCCESS;
}

