#pragma once
#include <bit>

namespace math
{
  class buffer
  {
  public:
    auto take();
    void free(unsigned const i);

  private:
    unsigned long long occupy = 0;
  };

  inline auto buffer::take()
  {
    auto index = std::countr_zero(~occupy);
    occupy |= (1ull << index);
    return index;
  }

  inline void buffer::free(unsigned const index)
  {
    occupy &= ~(1ull << index);
  }
}