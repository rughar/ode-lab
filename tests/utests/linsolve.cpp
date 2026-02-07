#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ling.hpp>

#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

// Deterministic sequence for test data.
// State is reduced modulo 1 and mapped to [-1,1) via a single-rounding FMA.
struct QuasiRandom
{
  double x = 0.0;
  static constexpr double phi = 0.6180339887498948482; // (sqrt(5)-1)/2 as a literal

  double next()
  {
    x += phi;
    x -= std::floor(x);            // x in [0,1)
    return std::fma(2.0, x, -1.0); // in [-1,1)
  }
};

static double residual_rel_inf(const std::vector<double> &A,
                               const std::vector<double> &y,
                               const std::vector<double> &b,
                               size_t n)
{
  std::vector<double> Ay(n, 0.0);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      Ay[i] += A[i * n + j] * y[j];

  double r = 0.0;
  for (size_t i = 0; i < n; ++i)
    r = std::max(std::abs(Ay[i] - b[i]), r);

  return r;
}

template <int D>
static void check_dim(double tol)
{
  static_assert(D >= 1 && D <= 10);
  constexpr size_t n = static_cast<size_t>(D);

  QuasiRandom rng;

  std::vector<double> M(D * D, 0.0), b(D, 0.0);

  for (size_t i = 0; i < D * D; ++i)
    M[i] = 0.2 * rng.next();

  for (size_t i = 0; i < D; ++i)
    b[i] = rng.next();

  for (size_t d = 0; d < D; ++d)
    M[d * D + d] += 1.0;

  // Prefix subproblem for this D.
  std::vector<double> A(M.begin(), M.begin() + D * D);
  std::vector<double> y(b.begin(), b.begin() + D);

  // In-place solve.
  math::solve_opt<double, D>(A.data(), y.data());

  const double rel = residual_rel_inf(M, y, b, n);

  INFO("D=" << D << " relative residual=" << rel);
  REQUIRE(rel < tol);
}

TEST_CASE("linsolve: solve_opt works on prefix subproblems", "[linsolve]")
{
  const double tol_rel = 1e-15;

  check_dim<1>(tol_rel);
  check_dim<2>(tol_rel);
  check_dim<3>(tol_rel);
  check_dim<4>(tol_rel);
  check_dim<5>(tol_rel);
  check_dim<6>(tol_rel);
  check_dim<7>(tol_rel);
  check_dim<8>(tol_rel);
  check_dim<9>(tol_rel);
  check_dim<10>(tol_rel);
}
