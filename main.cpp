#include <iostream>
#include <fstream>
#include <vector>

#include "include/predefined.h"

unsigned int read1(std::istream & in) {
  unsigned char a, b;
  a = static_cast<unsigned char>(in.get());
  b = static_cast<unsigned char>(in.get());
  return a | (b << 8);
}

std::vector<std::vector<std::size_t>> read_amalfi(std::istream & in) {
  std::size_t n = read1(in);
  std::vector<std::vector<std::size_t>> g(n);
  for(std::size_t i=0; i<n; ++i) {
    std::size_t cnt = read1(in);
    for(std::size_t j=0; j<cnt; ++j) {
      std::size_t k = read1(in);
      g[i].push_back(k);
    }
  }
  return g;
}

int main(int argc, char * argv[]) {
  char const * g_filename = argv[1];
  char const * h_filename = argv[2];

  std::ifstream in{g_filename,std::ios::in|std::ios::binary};
  auto g = read_amalfi(in);
  in.close();
  in.open(h_filename,std::ios::in|std::ios::binary);
  auto h = read_amalfi(in);
  in.close();
  
  int count = 0;

  ri_ind(
      g,
      h,
      [&count](auto const & S){++count; return true;},
      [](auto x, auto y){return true;},
      [](auto x0, auto x1, auto y0, auto y1){return true;});

  std::cout << count << std::endl;
}
