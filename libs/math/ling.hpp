// =============================================================================
//  FILE: ling.h  -  lightweight naive linear algebra utilities
// =============================================================================

#pragma once

#include <cmath>
#include <cstddef>

namespace math
{
  template <class U1, class U2>
  inline auto dot_product(const size_t n, const U1 a[], const U2 b[])
  {
    using U = std::remove_cvref_t<decltype(a[0] * b[0])>;
    U tmp = U(0);
    for (size_t i = 0; i < n; ++i)
      tmp += a[i] * b[i];
    return tmp;
  }

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

  template <class U>
  inline void lu_naive(const size_t n, U A[])
  {
    for (size_t i = 0; i < n; i++)
    {
      const size_t row_i = n * i;
      A[row_i + i] = 1 / A[row_i + i];
      const U a_ii = A[row_i + i];
      for (size_t j = i + 1; j < n; j++)
      {
        const size_t row_j = n * j;
        A[row_j + i] *= a_ii;
        const U a_ji = A[row_j + i];
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

  template <class U>
  inline void fb_naive(const size_t n, const U A[], U v[])
  {
    for (size_t i = 1; i < n; i++)
      v[i] -= dot_product(i, A + n * i, v);

    for (size_t i = n; i--;)
    {
      const size_t diag_i = (n + 1) * i;
      v[i] = A[diag_i] * (v[i] - dot_product(n - i - 1, A + diag_i + 1, v + i + 1));
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

  template <class U>
  inline auto spectral_radius_estimate(const size_t n, const U A[])
  {
    U tr1 = 0;
    U tr2 = 0;

    for (size_t i = 0; i < n; i++)
    {
      const size_t row_i = n * i;
      tr1 += A[row_i + i];
      for (size_t j = 0; j < n; j++)
        tr2 += A[row_i + j] * A[n * j + i];
    }

    const U det2 = 2 * tr2 - tr1 * tr1;

    if (det2 < 0)
      return std::sqrt(std::abs(tr2 / 2));

    return (std::abs(tr1) + std::sqrt(det2)) / 2;
  }

  template <class U, class Upack>
  inline void remove_tangent_components(const size_t n, const int p, U x[], const Upack &u)
  {
    for (int j = 0; j < p; j++)
    {
      U num = 0;

      for (size_t i = 0; i < n; ++i)
        num += u[j][i] * x[i];
      for (size_t i = 0; i < n; ++i)
        x[i] -= num * u[j][i];
    }
  }

  // ============================================================
  // solve_opt<U,n>(A,b)
  // - In-place solve of A*x = b with no pivoting.
  // - Overwrites A and b (b becomes x).
  // ============================================================

  template <class U, std::size_t n>
  inline void solve_opt(U A[n * n], U b[n])
  {
    if constexpr (n == 1)
    {
      // b0 = b0 / A00
      b[0] /= A[0];
    }
    else if constexpr (n == 2)
    {
      // Explicit 2x2 solve
      const U a00 = A[0], a01 = A[1];
      const U a10 = A[2], a11 = A[3];

      const U det = a00 * a11 - a01 * a10;
      const U invdet = U(1) / det;

      const U x0 = (a11 * b[0] - a01 * b[1]) * invdet;
      const U x1 = (-a10 * b[0] + a00 * b[1]) * invdet;

      b[0] = x0;
      b[1] = x1;
    }
    else if constexpr (n == 3)
    {
      // 3x3: single-shot elimination, store inv diag, back-sub with multiplies.
      // Layout: A[0..8] row-major.
      // k=0
      const U inv00 = U(1) / A[0];
      A[0] = inv00;

      const U f10 = A[3] * inv00;
      const U f20 = A[6] * inv00;

      // Update rows 1..2, cols 1..2 using row0 (A[1],A[2])
      A[4] -= f10 * A[1];
      A[5] -= f10 * A[2];

      A[7] -= f20 * A[1];
      A[8] -= f20 * A[2];

      b[1] -= f10 * b[0];
      b[2] -= f20 * b[0];

      // k=1
      const U inv11 = U(1) / A[4];
      A[4] = inv11;

      const U f21 = A[7] * inv11;

      A[8] -= f21 * A[5];
      b[2] -= f21 * b[1];

      // k=2
      const U inv22 = U(1) / A[8];
      A[8] = inv22;

      // back-sub
      b[2] *= A[8];
      b[1] -= A[5] * b[2];
      b[1] *= A[4];
      b[0] -= A[1] * b[1] + A[2] * b[2];
      b[0] *= A[0];
    }
    else if constexpr (n == 4)
    {
      // 4x4: as corrected one-shot with inv diag
      // k=0
      const U inv00 = U(1) / A[0];
      A[0] = inv00;

      const U f10 = A[4] * inv00;
      const U f20 = A[8] * inv00;
      const U f30 = A[12] * inv00;

      A[5] -= f10 * A[1];
      A[6] -= f10 * A[2];
      A[7] -= f10 * A[3];
      A[9] -= f20 * A[1];
      A[10] -= f20 * A[2];
      A[11] -= f20 * A[3];
      A[13] -= f30 * A[1];
      A[14] -= f30 * A[2];
      A[15] -= f30 * A[3];

      b[1] -= f10 * b[0];
      b[2] -= f20 * b[0];
      b[3] -= f30 * b[0];

      // k=1
      const U inv11 = U(1) / A[5];
      A[5] = inv11;

      const U f21 = A[9] * inv11;
      const U f31 = A[13] * inv11;

      A[10] -= f21 * A[6];
      A[11] -= f21 * A[7];
      A[14] -= f31 * A[6];
      A[15] -= f31 * A[7];

      b[2] -= f21 * b[1];
      b[3] -= f31 * b[1];

      // k=2
      const U inv22 = U(1) / A[10];
      A[10] = inv22;

      const U f32 = A[14] * inv22;

      A[15] -= f32 * A[11];
      b[3] -= f32 * b[2];

      // k=3
      const U inv33 = U(1) / A[15];
      A[15] = inv33;

      // back-sub
      b[3] *= A[15];

      b[2] -= A[11] * b[3];
      b[2] *= A[10];

      b[1] -= A[6] * b[2] + A[7] * b[3];
      b[1] *= A[5];

      b[0] -= A[1] * b[1] + A[2] * b[2] + A[3] * b[3];
      b[0] *= A[0];
    }
    else if constexpr (n == 5)
    {
      // 5x5: one-shot elimination, inv diag, no pivot.
      // Indices: row-major A[r*5 + c]

      // k=0
      const U inv00 = U(1) / A[0];
      A[0] = inv00;

      const U f10 = A[5] * inv00;
      const U f20 = A[10] * inv00;
      const U f30 = A[15] * inv00;
      const U f40 = A[20] * inv00;

      // update rows 1..4, cols 1..4 using row0 (A[1..4])
      A[6] -= f10 * A[1];
      A[7] -= f10 * A[2];
      A[8] -= f10 * A[3];
      A[9] -= f10 * A[4];
      A[11] -= f20 * A[1];
      A[12] -= f20 * A[2];
      A[13] -= f20 * A[3];
      A[14] -= f20 * A[4];
      A[16] -= f30 * A[1];
      A[17] -= f30 * A[2];
      A[18] -= f30 * A[3];
      A[19] -= f30 * A[4];
      A[21] -= f40 * A[1];
      A[22] -= f40 * A[2];
      A[23] -= f40 * A[3];
      A[24] -= f40 * A[4];

      b[1] -= f10 * b[0];
      b[2] -= f20 * b[0];
      b[3] -= f30 * b[0];
      b[4] -= f40 * b[0];

      // k=1
      const U inv11 = U(1) / A[6];
      A[6] = inv11;

      const U f21 = A[11] * inv11;
      const U f31 = A[16] * inv11;
      const U f41 = A[21] * inv11;

      // update rows 2..4, cols 2..4 using row1 (A[7..9])
      A[12] -= f21 * A[7];
      A[13] -= f21 * A[8];
      A[14] -= f21 * A[9];
      A[17] -= f31 * A[7];
      A[18] -= f31 * A[8];
      A[19] -= f31 * A[9];
      A[22] -= f41 * A[7];
      A[23] -= f41 * A[8];
      A[24] -= f41 * A[9];

      b[2] -= f21 * b[1];
      b[3] -= f31 * b[1];
      b[4] -= f41 * b[1];

      // k=2
      const U inv22 = U(1) / A[12];
      A[12] = inv22;

      const U f32 = A[17] * inv22;
      const U f42 = A[22] * inv22;

      // update rows 3..4, cols 3..4 using row2 (A[13..14])
      A[18] -= f32 * A[13];
      A[19] -= f32 * A[14];
      A[23] -= f42 * A[13];
      A[24] -= f42 * A[14];

      b[3] -= f32 * b[2];
      b[4] -= f42 * b[2];

      // k=3
      const U inv33 = U(1) / A[18];
      A[18] = inv33;

      const U f43 = A[23] * inv33;

      // update row 4, col 4 using row3 (A[19])
      A[24] -= f43 * A[19];
      b[4] -= f43 * b[3];

      // k=4
      const U inv44 = U(1) / A[24];
      A[24] = inv44;

      // back-sub
      b[4] *= A[24];

      b[3] -= A[19] * b[4];
      b[3] *= A[18];

      b[2] -= A[13] * b[3] + A[14] * b[4];
      b[2] *= A[12];

      b[1] -= A[7] * b[2] + A[8] * b[3] + A[9] * b[4];
      b[1] *= A[6];

      b[0] -= A[1] * b[1] + A[2] * b[2] + A[3] * b[3] + A[4] * b[4];
      b[0] *= A[0];
    }
    else if constexpr (n == 6)
    {
      // 6x6: corrected one-shot with inv diag stored
      // k=0
      const U inv00 = U(1) / A[0];
      A[0] = inv00;

      const U f10 = A[6] * inv00;
      const U f20 = A[12] * inv00;
      const U f30 = A[18] * inv00;
      const U f40 = A[24] * inv00;
      const U f50 = A[30] * inv00;

      A[7] -= f10 * A[1];
      A[8] -= f10 * A[2];
      A[9] -= f10 * A[3];
      A[10] -= f10 * A[4];
      A[11] -= f10 * A[5];
      A[13] -= f20 * A[1];
      A[14] -= f20 * A[2];
      A[15] -= f20 * A[3];
      A[16] -= f20 * A[4];
      A[17] -= f20 * A[5];
      A[19] -= f30 * A[1];
      A[20] -= f30 * A[2];
      A[21] -= f30 * A[3];
      A[22] -= f30 * A[4];
      A[23] -= f30 * A[5];
      A[25] -= f40 * A[1];
      A[26] -= f40 * A[2];
      A[27] -= f40 * A[3];
      A[28] -= f40 * A[4];
      A[29] -= f40 * A[5];
      A[31] -= f50 * A[1];
      A[32] -= f50 * A[2];
      A[33] -= f50 * A[3];
      A[34] -= f50 * A[4];
      A[35] -= f50 * A[5];

      b[1] -= f10 * b[0];
      b[2] -= f20 * b[0];
      b[3] -= f30 * b[0];
      b[4] -= f40 * b[0];
      b[5] -= f50 * b[0];

      // k=1
      const U inv11 = U(1) / A[7];
      A[7] = inv11;

      const U f21 = A[13] * inv11;
      const U f31 = A[19] * inv11;
      const U f41 = A[25] * inv11;
      const U f51 = A[31] * inv11;

      A[14] -= f21 * A[8];
      A[15] -= f21 * A[9];
      A[16] -= f21 * A[10];
      A[17] -= f21 * A[11];
      A[20] -= f31 * A[8];
      A[21] -= f31 * A[9];
      A[22] -= f31 * A[10];
      A[23] -= f31 * A[11];
      A[26] -= f41 * A[8];
      A[27] -= f41 * A[9];
      A[28] -= f41 * A[10];
      A[29] -= f41 * A[11];
      A[32] -= f51 * A[8];
      A[33] -= f51 * A[9];
      A[34] -= f51 * A[10];
      A[35] -= f51 * A[11];

      b[2] -= f21 * b[1];
      b[3] -= f31 * b[1];
      b[4] -= f41 * b[1];
      b[5] -= f51 * b[1];

      // k=2
      const U inv22 = U(1) / A[14];
      A[14] = inv22;

      const U f32 = A[20] * inv22;
      const U f42 = A[26] * inv22;
      const U f52 = A[32] * inv22;

      A[21] -= f32 * A[15];
      A[22] -= f32 * A[16];
      A[23] -= f32 * A[17];
      A[27] -= f42 * A[15];
      A[28] -= f42 * A[16];
      A[29] -= f42 * A[17];
      A[33] -= f52 * A[15];
      A[34] -= f52 * A[16];
      A[35] -= f52 * A[17];

      b[3] -= f32 * b[2];
      b[4] -= f42 * b[2];
      b[5] -= f52 * b[2];

      // k=3
      const U inv33 = U(1) / A[21];
      A[21] = inv33;

      const U f43 = A[27] * inv33;
      const U f53 = A[33] * inv33;

      A[28] -= f43 * A[22];
      A[29] -= f43 * A[23];
      A[34] -= f53 * A[22];
      A[35] -= f53 * A[23];

      b[4] -= f43 * b[3];
      b[5] -= f53 * b[3];

      // k=4
      const U inv44 = U(1) / A[28];
      A[28] = inv44;

      const U f54 = A[34] * inv44;

      A[35] -= f54 * A[29];
      b[5] -= f54 * b[4];

      // k=5
      const U inv55 = U(1) / A[35];
      A[35] = inv55;

      // back-sub
      b[5] *= A[35];

      b[4] -= A[29] * b[5];
      b[4] *= A[28];

      b[3] -= A[22] * b[4] + A[23] * b[5];
      b[3] *= A[21];

      b[2] -= A[15] * b[3] + A[16] * b[4] + A[17] * b[5];
      b[2] *= A[14];

      b[1] -= A[8] * b[2] + A[9] * b[3] + A[10] * b[4] + A[11] * b[5];
      b[1] *= A[7];

      b[0] -= A[1] * b[1] + A[2] * b[2] + A[3] * b[3] + A[4] * b[4] + A[5] * b[5];
      b[0] *= A[0];
    }
    else
    {
      lu_naive(n, A);
      fb_naive(n, A, b);
    }
  }

}