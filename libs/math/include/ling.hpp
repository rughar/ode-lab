// =============================================================================
//  FILE: ling.h  -  lightweight naive LU solver
// =============================================================================
//  Purpose :
//    * lu_naive  - factorises a diagonally dominant matrix (1+A) into L*U.
//    * fb_naive  - solves (1+A)*x = b in place on v by forward + backward
//      substitution using the factors produced by lu_naive.
//
//  Limitations :
//    - No pivoting (numeric stability relies on diagonal dominance).
//    - Works on user supplied containers that support A[n * i + j] style access
//      and are mutable for the LU stage.
// =============================================================================

#pragma once

namespace math
{

  // ---------------------------------------------------------------------------
  //  lu_naive
  // ---------------------------------------------------------------------------
  //  In‑place LU factorisation of an  n*n  matrix  (1+A)  without pivoting.
  //  After the call:
  //      • U  is stored on and above the main diagonal.
  //      • L  (without the unit diagonal) is stored below the diagonal.
  //  PRECONDITION:  A  must be diagonally dominant so that no element on the
  //  diagonal becomes zero.
  // ---------------------------------------------------------------------------

  template <class Um>
  void lu_naive(const size_t n, Um &A)
  {
    for (size_t i = 0; i < n; i++)
    {
      A[n * i + i] = 1 / A[n * i + i];
      auto aii = A[n * i + i];
      for (size_t j = i + 1; j < n; j++)
      {
        A[n * j + i] *= aii;
        auto aji = A[n * j + i];
        for (size_t k = i + 1; k < n; k++)
          A[n * j + k] -= aji * A[n * i + k];
      }
    }
  }

  // ---------------------------------------------------------------------------
  //  fb_naive
  // ---------------------------------------------------------------------------
  //  Forward + backward substitution for  (1+A).x = v .
  //  The matrix  (1+A)  must already be factorised by  lu_naive.
  //  The solution overwrites  v.
  // ---------------------------------------------------------------------------

  template <class Um, class Uv>
  void fb_naive(const size_t n, const Um &A, Uv &v)
  {
    for (size_t i = 1; i < n; i++)
      for (size_t j = 0; j < i; j++)
        v[i] -= A[n * i + j] * v[j];
    for (size_t i = n; i--;)
    {
      for (size_t j = i + 1; j < n; j++)
        v[i] -= A[n * i + j] * v[j];
      v[i] *= A[n * i + i];
    }
  }
}