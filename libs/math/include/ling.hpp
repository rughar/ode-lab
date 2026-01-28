// =============================================================================
//  FILE: ling.h  -  lightweight naive linear algebra utilities
// =============================================================================
//  Interface :
//    - Works on user-supplied containers that support A[n*i + j] indexing
//      via operator[].
// =============================================================================

#pragma once

#include <cmath>

namespace math
{

  // ---------------------------------------------------------------------------
  //  lu_naive
  // ---------------------------------------------------------------------------
  //  In‑place LU factorisation of an n*n matrix A without pivoting.
  //  After the call:
  //      • U  is stored on and above the main diagonal.
  //      • L  (without the unit diagonal) is stored below the diagonal.
  //  PRECONDITION: A must be diagonally dominant so that no element on the
  //  diagonal becomes zero.
  // ---------------------------------------------------------------------------

  template <class Um>
  inline void lu_naive(const size_t n, Um &A)
  {
    for (size_t i = 0; i < n; i++)
    {
      const size_t row_i = n * i;
      A[row_i + i] = 1 / A[row_i + i];
      auto a_ii = A[row_i + i];
      for (size_t j = i + 1; j < n; j++)
      {
        const size_t row_j = n * j;
        A[row_j + i] *= a_ii;
        auto a_ji = A[row_j + i];
        for (size_t k = i + 1; k < n; k++)
          A[row_j + k] -= a_ji * A[row_i + k];
      }
    }
  }

  // ---------------------------------------------------------------------------
  //  fb_naive
  // ---------------------------------------------------------------------------
  //  Forward + backward substitution for A.x = v
  //  The matrix A must already be factorised by  lu_naive.
  //  The solution overwrites v
  // ---------------------------------------------------------------------------

  template <class Um, class Uv>
  inline void fb_naive(const size_t n, const Um &A, Uv &v)
  {
    for (size_t i = 1; i < n; i++)
    {
      const size_t row_i = n * i;
      for (size_t j = 0; j < i; j++)
        v[i] -= A[row_i + j] * v[j];
    }

    for (size_t i = n; i--;)
    {
      const size_t row_i = n * i;
      for (size_t j = i + 1; j < n; j++)
        v[i] -= A[row_i + j] * v[j];
      v[i] *= A[row_i + i];
    }
  }

  // ---------------------------------------------------------------------------
  //  spectral_radius_estimate
  // ---------------------------------------------------------------------------
  //  Estimates the magnitude of the dominant eigenvalue (spectral radius) from
  //  two quadratic invariants of the matrix:
  //
  //      • tr1 = trace(A)
  //      • tr2 = trace(A.A)
  //
  //  return value: max(|λ1|, |λ2|)
  //
  //  For n == 2 the eigenvalues λ1, λ2 satisfy:
  //
  //      λ1 + λ2 = tr1
  //      λ1^2 + λ2^2 = tr2
  //
  //  (i.e. the spectral radius for a real 2×2 spectrum).
  //  Note: For n > 2 this is only a heuristic but surprisingly good estimate
  //  where we expect that after matrix multiply A := A.A smaller eigenvalues
  //  will be less dominant.
  // ---------------------------------------------------------------------------

  template <class Um>
  inline auto spectral_radius_estimate(const size_t n, const Um &A)
  {
    using U = std::remove_cvref_t<decltype(A[0])>;

    U tr1 = 0;
    U tr2 = 0;

    for (size_t i = 0; i < n; i++)
    {
      const size_t row_i = n * i;
      tr1 += A[row_i + i];
      for (size_t j = 0; j < n; j++)
        tr2 += A[row_i + j] * A[n * j + i];
    }

    U det2 = 2 * tr2 - tr1 * tr1;

    if (det2 < 0)
      return std::sqrt(std::abs(tr2 / 2));

    return (std::abs(tr1) + std::sqrt(det2)) / 2;
  }

  
  template<class V1, class V2>
  inline auto dot_product(const size_t n, const V1 &a, const V2 &b)
  {
    using U = std::remove_cvref_t<decltype(a[0] * b[0])>;
    U tmp = 0;
    for(size_t i = 0; i < n; ++i)
      tmp += a[i] * b[i];
    return tmp;
  }


  template <class V, class Vpack>
  inline void remove_tangent_components(const size_t n, const int p, V &x, const Vpack &u)
  {
    for (int j = 0; j < p; j++)
    {
      using U = std::remove_cvref_t<decltype(x[0])>;
      U num = 0;
      
      for (size_t i = 0; i < n; ++i)
        num += u[j][i] * x[i];
      for (size_t i = 0; i < n; ++i)
        x[i] -= num * u[j][i];
    }
  }

}