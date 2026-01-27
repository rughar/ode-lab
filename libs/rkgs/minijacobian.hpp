#pragma once

#include <vector>

namespace rkgs
{
  template<class U, int level>
  class mini_jacobian 
  {
    public:
      template<class F, class V>
      void evaluate(F& f, const V& x, const U radius);
    private:
      size_t n;
      std::vector<U> u[level], v[level];

  };
  
  template <class U, int level>
  template <class F, class V>
  inline void mini_jacobian<U, level>::evaluate(F &f, const V &x, const U radius)
  {
  }
}