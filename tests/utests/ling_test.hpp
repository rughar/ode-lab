#pragma once
#include <utest_frame.hpp>
#include <ling.hpp>
#include <string>

void test_dot_product()
{
  double a[] = {1.0, 2.0, 3.0};
  double b[] = {4.0, 5.0, 6.0};
  double expected = 1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0;
  utest::compare_numeric("wrong dot product result", expected, math::dot_product(3, a, b));
}

void test_spectral_radius_estimate()
{
  const double mat[] = {
      1.1,
      2.0,
      0.0,
      0.8,
  };

  double expected = 1.1;
  double result = math::spectral_radius_estimate(2, mat);
  utest::compare_numeric("wrong spectral radius on 2*2 matrix", expected, result);

  const double mat2[] = {0.0, 3.0, -3.0, 0.0};
  expected = 3.0;
  result = math::spectral_radius_estimate(2, mat2);
  utest::compare_numeric("wrong spectral radius on 2*2 antisymmetric matrix", expected, result);
}

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

template <size_t n>
void subtest_solve_opt(utest::error_accumulator &acc)
{
  double A[n * n], B[n * n];
  double x[n], y[n];

  QuasiRandom qr;

  for (size_t i = 0; i < n * n; ++i)
    A[i] = B[i] = 0.1 * qr.next();
  for (size_t i = 0; i < n; ++i)
  {
    x[i] = y[i] = qr.next();
    A[n * i + i] += 1.0;
    B[n * i + i] += 1.0;
  }

  math::solve_opt<n>(A, x);

  for (size_t i = 0; i < n; ++i)
  {
    double sum = math::dot_product(n, B + n * i, x);
    try {
      utest::compare_numeric("wrong solve_opt<" + std::to_string(n) + ">", y[i], sum, 2e-16 * std::sqrt(n));
    }
    catch (const std::exception &e)
    {
      acc.add(e.what());
    }
  }
}

void test_solve_opt()
{
  utest::error_accumulator acc;

  subtest_solve_opt<1>(acc);
  subtest_solve_opt<2>(acc);
  subtest_solve_opt<3>(acc);
  subtest_solve_opt<4>(acc);
  subtest_solve_opt<5>(acc);
  subtest_solve_opt<6>(acc);
  subtest_solve_opt<7>(acc);
  subtest_solve_opt<20>(acc);

  acc.throw_if_any();
}

void test_remove_tangent_components()
{
  double u[2][3] = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
  };
  double x[3] = {1.0, 2.0, 3.0};
  math::remove_tangent_components(3, 2, x, u);
  utest::compare_numeric("wrong remove_tangent_component 0", 0.0, x[0]);
  utest::compare_numeric("wrong remove_tangent_component 1", 0.0, x[1]);
  utest::compare_numeric("wrong remove_tangent_component 2", 3.0, x[2]);
}