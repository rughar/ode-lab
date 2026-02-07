#pragma once
#include <vector>
#include <minijacobian.hpp>

namespace rkgl
{
  template <class U, int order>
  inline const U A[order + 1][order] = {};

  template <class U>
  inline const U A<U, 1>[2][1] = {
      {0.5},
      {1.0}};

#define rsqrt(a) (1 / std::sqrt(a))

  template <class U>
  inline const U A<U, 2>[3][2] = {
      {0.25, 0.25 - 0.5 * rsqrt(3.0)},
      {0.25 + 0.5 * rsqrt(3.0), 0.25},
      {0.5, 0.5}};

  template <class U>
  inline const U A<U, 3>[4][3] = {
      {5.0 / 36.0, 2.0 / 9.0 - rsqrt(15.0), 5.0 / 36.0 - 0.5 * rsqrt(15.0)},
      {5.0 / 36.0 + 0.625 * rsqrt(15.0), 2.0 / 9.0, 5.0 / 36.0 - 0.625 * rsqrt(15.0)},
      {5.0 / 36.0 + 0.5 * rsqrt(15.0), 2.0 / 9.0 + rsqrt(15.0), 5.0 / 36.0},
      {5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0}};

  template <class U, int order>
  class rkgl
  {
  public:
    template <class F>
    void step(F &f, mini_jacobian<U> jac, U *x, const U h, const U tol);
    void set(size_t size);

  private:
    size_t n;
    std::vector<U> k[order], z, y, k_tmp;
  };

  template <class U, int order>
  template <class F>
  inline void rkgl<U, order>::step(F &f, mini_jacobian<U> jac, U *x, const U h, const U tol)
  {
    for (int s = 0; s < 100; s++)
    {
      for (int j = 0; j < order; j++)
      {
        for (size_t i = 0; i < n; i++)
        {
          U tmp = U(0);
          for (int m = 0; m < order; m++)
            tmp += k[m][i] * A<U, order>[j][m];
          z[i] = x[i] + h * tmp;
        }

        k_tmp = k[j];
        f(z.data(), k[j].data());

        for (size_t i = 0; i < n; i++)
          k_tmp[i] = h * (k[j][i] - k_tmp[i]) ;
        jac.aply(k_tmp.data(), k_tmp.data());

        for (int m = 0; m < order; m++)
          for (size_t i = 0; i < n; i++)
            k[m][i] += k_tmp[i] * A<U, order>[m][j];
      }



    }
  };

  template <class U, int order>
  inline void rkgl<U, order>::set(size_t size)
  {
    n = size;
    for (int j = 0; j < order; j++)
      k[j].resize(n, U(0));
    y.resize(n, U(0));
    z.resize(n, U(0));
  }
}
