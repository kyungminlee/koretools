#include <iostream>

#include <kore/kore>
#include <kore/utility/tuple_stream.h>
#include <kore/bitbox/bitshifter.h>


int main(int argc, char** argv)
{
  using namespace kore::bitbox;

  int64_t x=2,y=3,z=4;
  BitShifter<3, uint64_t, int64_t> bs(x,y,z);
  std::cout << bs.size() << std::endl;
  std::cout << std::bitset<32>(bs.mask_all()) << std::endl;
  std::cout << bs.shape() << std::endl;
  std::cout << bs.stride() << std::endl;
  std::cout << std::bitset<32>(bs.mask()[0]) << std::endl;
  std::cout << std::bitset<32>(bs.mask()[1]) << std::endl;
  std::cout << std::bitset<32>(bs.mask()[2]) << std::endl;

  return 0;
}

