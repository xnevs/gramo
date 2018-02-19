#include <iostream>
#include <fstream>
#include <vector>

#include "include/read_amalfi.h"
#include "include/simple_adjacency_list.h"
#include "include/predefined.h"

int main(int argc, char * argv[]) {
  char const * g_filename = argv[1];
  char const * h_filename = argv[2];

  std::ifstream in{g_filename,std::ios::in|std::ios::binary};
  auto g = read_amalfi<simple_adjacency_list<uint16_t>>(in);
  in.close();
  in.open(h_filename,std::ios::in|std::ios::binary);
  auto h = read_amalfi<simple_adjacency_list<uint16_t>>(in);
  in.close();
  
  
  int count = 0;
  
  //ullmann_ind_RDEG_CNC(
  //ri_ind(
  ri_ind(
      g,
      h,
      [&count](auto const & S) {++count; return true;},
      [](auto x, auto y) {return true;},
      [](auto x0, auto x1, auto y0, auto y1) {return true;});

  std::cout << count << std::endl;
}
