#pragma once
#include <utest_frame.hpp>
#include <ling.hpp>
#include <string>

void test_dot_product()
{
  double a[] = {1.0, 2.0, 3.0};
  double b[] = {4.0, 5.0, 6.0};
  double expected = 1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0;
  double result = math::dot_product(3, a, b);
  if (result != expected)
    throw std::runtime_error("wrong dot product result");
}

void test_spectral_radius_estimate()
{
  double mat[] = {
      1.0,
      2.0,
      0.0,
      1.0,
  };
  double expected = 1.0;
  double result = math::spectral_radius_estimate(2, mat);
  if (result != expected)
    throw std::runtime_error("wrong spectral radius on 2*2 matrix");
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
void subtest_solve_opt()
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
    if (std::abs(sum - y[i]) > 1e-15)
      throw std::runtime_error("wrong solve_opt<" + std::to_string(n) + ">");
  }
}

void test_solve_opt()
{
  subtest_solve_opt<1>();
  subtest_solve_opt<2>();
  subtest_solve_opt<3>();
  subtest_solve_opt<4>();
  subtest_solve_opt<5>();
  subtest_solve_opt<6>();
  subtest_solve_opt<7>();
  subtest_solve_opt<8>();
}

void test_remove_tangent_components()
{
  double u[2][3] = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
  };
  double x[3] = {1.0, 2.0, 3.0};
  math::remove_tangent_components(3, 2, x, u);
  if (std::abs(x[0]) > 1e-15 || std::abs(x[1]) > 1e-15 || std::abs(x[2] - 3.0) > 1e-15)
    throw std::runtime_error("wrong remove_tangent_components");
}