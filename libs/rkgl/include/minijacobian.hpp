#pragma once

#include <vector>
#include <ling.hpp>
#include <buffer.hpp>

namespace rkgl
{
  template <class U>
  class mini_jacobian
  {
  public:
    template <class F>
    int evaluate(F &f, const U *x, const U h);

    template <class F, class CV>
    int evaluate(F &f, const U *x, const U h, CV &&covector);

    void aply(const U *x, U *y) const;
    U spectral_radius_estimate() const;

    void set(size_t size);

  protected:
    size_t n;

  private:
    std::vector<U> u[2], v[2], y0, buf[2];
    U udotv[4];

    static inline constexpr auto euclid_covector = [](size_t n, const U *x, U *cx)
    {
      for (size_t i = 0; i < n; ++i)
        cx[i] = x[i];
    };
  };

  template <class U>
  template <class F>
  inline int mini_jacobian<U>::evaluate(F &f, const U *x, const U h)
  {
    return evaluate(f, x, h, euclid_covector);
  }

  template <class U>
  template <class F, class CV>
  inline int mini_jacobian<U>::evaluate(F &f, const U *x, const U h, CV &&covector)
  {
    f(x, y0.data());

    auto &dx = buf[0];
    auto &x_push = buf[1];
    
    for (size_t i = 0; i < n; ++i)
      dx[i] = y0[i] * h;

    for (int k : {0, 1})
    {
      auto &pu = u[k];
      auto &pv = v[k];

      for (size_t i = 0; i < n; ++i)
        x_push[i] = x[i] + dx[i];
      f(x_push.data(), pv.data());

      for (size_t i = 0; i < n; ++i)
        x_push[i] -= 2 * dx[i];
      f(x_push.data(), pu.data());

      auto &pw = x_push;

      for (size_t i = 0; i < n; ++i)
      {
        pw[i] = (pv[i] + pu[i]) / 2;
        pv[i] = (pv[i] - pu[i]) / 2;
      }

      covector(n, dx.data(), pu.data());
      U denom = math::dot_product(n, pu.data(), dx.data());

      if (denom == U(0))
        return k;

      denom = 1 / denom;
      U alpha_v = udotv[3 * k] = denom * math::dot_product(n, pu.data(), pv.data());

      if (k == 0)
      {
        U alpha_w = denom * math::dot_product(n, pu.data(), pw.data());

        auto &pa = u[1];
        auto &pb = v[1];
        auto &pca = buf[0];
        auto &pcb = buf[1];

        for (size_t i = 0; i < n; ++i)
        {
          pa[i] = pw[i] - alpha_w * dx[i];
          pb[i] = pv[i] - alpha_v * dx[i];
        }

        covector(n, pa.data(), pca.data());
        covector(n, pb.data(), pcb.data());

        if (math::dot_product(n, pa.data(), pca.data()) > math::dot_product(n, pb.data(), pcb.data()))
          for (size_t i = 0; i < n; ++i)
            dx[i] = h * pa[i];
        else
          for (size_t i = 0; i < n; ++i)
            dx[i] = h * pb[i];
      }

      denom = std::sqrt(denom);
      for (size_t i = 0; i < n; ++i)
      {
        pu[i] *= denom;
        pv[i] *= denom;
      }
    }

    udotv[1] = math::dot_product(n, u[0].data(), v[1].data());
    udotv[2] = math::dot_product(n, u[1].data(), v[0].data());

    return 2;
  }

  template <class U>
  inline void mini_jacobian<U>::aply(const U x[], U y[]) const
  {
    U p[2] = {math::dot_product(n, u[0].data(), x), math::dot_product(n, u[1].data(), x)};
    for (size_t i = 0; i < n; ++i)
      y[i] = p[0] * v[0][i] + p[1] * v[1][i];
  }

  template <class U>
  inline U mini_jacobian<U>::spectral_radius_estimate() const
  {
    return math::spectral_radius_estimate(2, udotv);
  }

  template <class U>
  inline void mini_jacobian<U>::set(size_t size)
  {
    n = size;
    u[0].resize(n);
    u[1].resize(n);
    v[0].resize(n);
    v[1].resize(n);
    buf[0].resize(n);
    buf[1].resize(n);
    y0.resize(n);
  }
}