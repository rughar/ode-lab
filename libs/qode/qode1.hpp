#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <ling.hpp>

// =============================================================================
//  FILE: qode1.hpp  -  Quadratic ODE integrator with stepsize control
// =============================================================================
//
//  Purpose
//  -------
//  This file defines the class qode1_core<U>, which provides a symmetric
//  one-step integrator for systems of quadratic ordinary differential equations
//
//      \dot{x}_i = A_i
//                  + sum_j     B_{i,j} x_j
//                  + sum_{j,k} C_{i,j,k} x_j x_k
//
//  The system dimension is specified at construction time. The state vector
//  is stored in the public member `u` and is updated in place.
//
//
//  Integration scheme
//  ------------------
//  Each step performs an implicit symmetric update of the form
//
//      x_{n+1} = x_n + h * ( A + B . (x_n + x_{n+1})/2 + C . x_n . x_n{n+1} ),
//
//  The resulting linear system is solved using LU factorisation.
//
//
//  How to use
//  ----------
//  1) Derive a class from qode1_core<U>.
//  2) Pass the system dimension to the base-class constructor.
//  3) Override the method
//
//        void set_coef();
//
//     Inside set_coef(), assign coefficients using:
//
//        A_COEF(i)        = value;
//        B_COEF(i, j)     = value;
//        C_COEF(i, j, k)  = value;
//
//  4) Set the initial state in the public vector `u`.
//  5) Use
//
//        step(h)
//
//     for step with fixed stepsize h.
//
//
//  Symmetric stepsize control
//  --------------------------
//  The class provides stepsize control steps based on an estimate of
//  the Jacobian spectral radius:
//
//    * h_start = suggest_first_stepsize(h_max, mu)
//
//        Suggests an initial stepsize h_start such that
//            spectral_radius * h_start = mu,
//        limited by the maximal allowed stepsize h_max.
//
//    * step_adaptive(h, mu, low_bound, high_bound)
//
//        Adjusts h so that the dimensionless stability measure
//            spectral_radius * h = mu
//        and performs a step with the adjusted stepsize.
//        The parameters low_bound (0 < low_bound ≤ 1) and
//        high_bound (high_bound ≥ 1) limit the minimal and maximal
//        multiplicative change of h relative to its previous value.
//
//        The stepsize selection accounts for time symmetry, meaning that two
//        consecutive steps h_old and h_new satisfy the relation
//            h_new * h_old = h_mid^2,
//        where h_mid is defined by
//            spectral_radius * h_mid = mu.
//
//
//  Notes
//  -----
//  * The state vector `u` is updated in place.
//  * No pivoting is used in the linear solver; numerical stability relies
//    on moderate stepsizes and problem structure.
//  * The Jacobian spectral radius for stepsize is estimated heuristically and
//    is not guaranteed to be an upper or lower bound.
//
// =============================================================================

namespace qode
{
  template <class U>
  class qode1_core
  {
  public:
    std::vector<U> x;

    explicit qode1_core(const size_t size);

    virtual void set_coef() = 0;

    size_t dim() const;
    void step(const U h);
    void step_adaptive(U &h, const U mu, const U low_bound = U(0.3), const U high_bound = U(2.0));
    U suggest_first_stepsize(const U h_max, const U mu);

  protected:
    std::vector<U> mat, vec;

    struct ACoefProxy
    {
      qode1_core &self;
      size_t i;
      void operator=(U value);
    };

    struct BCoefProxy
    {
      qode1_core &self;
      size_t i, j;
      void operator=(U value);
    };

    struct CCoefProxy
    {
      qode1_core &self;
      size_t i, j, k;
      void operator=(U value);
    };

    ACoefProxy a_coef(const size_t i);
    BCoefProxy b_coef(const size_t i, const size_t j);
    CCoefProxy c_coef(const size_t i, const size_t j, const size_t k);

  private:
    size_t n;

    void prepare_step();
    void finish_step(const U h);
    U jacobian_spectral_radius();
  };

  // -------------------------------------------------------------------------
  //  qode1_core<U> implementation
  // -------------------------------------------------------------------------

  // -- public API ------------------------------------------------------------

  template <class U>
  inline qode1_core<U>::qode1_core(const size_t size) : n(size)
  {
    vec.resize(n, U(0));
    mat.resize(n * n, U(0));
  }

  template <class U>
  inline size_t qode1_core<U>::dim() const
  {
    return n;
  }

  template <class U>
  inline void qode1_core<U>::step(const U h)
  {
    prepare_step();
    finish_step(h);
  }

  template <class U>
  inline void qode1_core<U>::step_adaptive(U &h, const U mu, const U low_bound, const U high_bound)
  {
    prepare_step();
    U omega = jacobian_spectral_radius();
    h *= std::max(low_bound, std::sqrt(mu / std::max(mu / (high_bound * high_bound), omega * h)));
    finish_step(h);
  }

  template <class U>
  inline U qode1_core<U>::suggest_first_stepsize(const U h_max, const U mu)
  {
    prepare_step();
    U omega = jacobian_spectral_radius();
    return mu / std::max(mu / h_max, omega);
  }

  // -- proxies ---------------------------------------------------------------

  template <class U>
  inline void qode1_core<U>::ACoefProxy::operator=(U value)
  {
    self.vec[i] += value;
  }

  template <class U>
  inline void qode1_core<U>::BCoefProxy::operator=(U value)
  {
    self.mat[self.n * i + j] += value;
    self.vec[i] += value * self.x[j] / 2;
  }

  template <class U>
  inline void qode1_core<U>::CCoefProxy::operator=(U value)
  {
    self.mat[self.n * i + j] += value * self.x[k];
    self.mat[self.n * i + k] += value * self.x[j];
  }

  template <class U>
  inline typename qode1_core<U>::ACoefProxy
  qode1_core<U>::a_coef(const size_t i)
  {
    return ACoefProxy{*this, i};
  }

  template <class U>
  inline typename qode1_core<U>::BCoefProxy
  qode1_core<U>::b_coef(const size_t i, const size_t j)
  {
    return BCoefProxy{*this, i, j};
  }

  template <class U>
  inline typename qode1_core<U>::CCoefProxy
  qode1_core<U>::c_coef(const size_t i, const size_t j, const size_t k)
  {
    return CCoefProxy{*this, i, j, k};
  }

  // -- private ---------------------------------------------------------------

  template <class U>
  inline U qode1_core<U>::jacobian_spectral_radius()
  {
    return math::spectral_radius_estimate(n, mat.data());
  }

  template <class U>
  inline void qode1_core<U>::prepare_step()
  {
    std::fill(vec.begin(), vec.end(), U(0));
    std::fill(mat.begin(), mat.end(), U(0));
    set_coef();
  }

  template <class U>
  inline void qode1_core<U>::finish_step(const U h)
  {
    for (size_t i = 0; i < n; i++)
    {
      const size_t row_i = n * i;
      for (size_t j = 0; j < n; j++)
        mat[row_i + j] *= -h / 2;

      mat[row_i + i] += 1;
      x[i] += h * vec[i];
    }

    math::lu_naive(n, mat.data());
    math::fb_naive(n, mat.data(), x.data());
  }
}