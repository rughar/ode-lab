#pragma once

#include <vector>
#include <ling.hpp>

namespace rkgs
{
  template <class U, int level>
  class mini_jacobian
  {
  public:

    template <class F, class V> void evaluate(F &f, const V &x, const U radius);
    template <class V1, class V2> void aply(const V1 &dx, V2 &dy) const;
    U spectral_radius_estimate();

  protected:

    size_t n;
    template <class V1, class V2, class V3> virtual void covector(const V1 &src, V2 &dst, const V3& x) const;

  private:

    std::vector<U> u[level], v[level], dx, y0, z;
    U dp[level * level];

  };

  template <class U, int level>
  template <class F, class V>
  inline void mini_jacobian<U, level>::evaluate(F &f, const V &x, const U r)
  {
    f(x, y0);

    for (size_t i = 0; i < n; ++i)
      dx[i] = y0[i] * r;

    for (int k = 0; k < level; ++k)
    {
      math::remove_tangent_components(n, k, dx, u);

      for (size_t i = 0; i < n; ++i)
        z[i] = x[i] + dx[i];

      f(z, v[k]);

      for (size_t i = 0; i < n; ++i)
        z[i] = x[i] - dx[i];

      f(z, u[k]);

      for (size_t i = 0; i < n; ++i)
        v[k][i] = (v[k][i] - u[k][i]) / 2;

      covector( dx, u[k] );
      U tmp = 1 / math::dot_product( n, u[k], dx );
    
      for (size_t i = 0; i < n; ++i)
        u[k][i] *= tmp;
    
      
    }
  }

  template <class U, int level>
  template <class V1, class V2>
  inline void mini_jacobian<U, level>::aply(const V1 &dx, V2 &dy) const
  {
    for(size_t i = 0; i < n; i++) {
      //dy[i]
    }
  }

  template <class U, int level>
  template <class V1, class V2, class V3>
  inline void mini_jacobian<U, level>::covector(const V1 &src, V2 &dst, const V3 &x) const
  {
    for(size_t i = 0; i < n; i++) {
      dst[i] = src[i];
  }

  template <class U, int level>
  inline U mini_jacobian<U, level>::spectral_radius_estimate()
  {
    if (level == 1)
      return dot_product(u[0], v[0]);

    for (int i = 0; i < level; ++i)
      for (int j = 0; j < level; ++j)
        dp[level * i + j] = dot_product(u[i], v[j]);

    return math::spectral_radius_estimate(level, dp);
  }
}