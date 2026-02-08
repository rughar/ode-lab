#pragma once
#include <utest_frame.hpp>
#include <ling.hpp>
#include <string>

template<class U>
constexpr U eps = std::numeric_limits<U>::epsilon();

void test_dot_product(utest::error_accumulator &ea)
{
  double a[] = {1.0, 2.0, 3.0};
  double b[] = {4.0, 5.0, 6.0};
  double expected = 1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0;
  ea << utest::compare_numeric("wrong dot product result", expected, math::dot_product(3, a, b));
}

void test_spectral_radius_estimate(utest::error_accumulator &ea)
{
  const double mat[] = {
      1.1,
      2.0,
      0.0,
      0.8,
  };

  double expected = 1.1;
  double result = math::spectral_radius_estimate(2, mat);
  ea << utest::compare_numeric("wrong spectral radius on 2*2 matrix", expected, result, eps<double>);

  const double mat2[] = {0.0, 3.0, -3.0, 0.0};
  expected = 3.0;
  result = math::spectral_radius_estimate(2, mat2);
  ea << utest::compare_numeric("wrong spectral radius on 2*2 antisymmetric matrix", expected, result, eps<double>);
}

struct QuasiRandom
{
  double x = 0.0;
  static constexpr double phi = 0.6180339887498948482; // (sqrt(5)-1)/2

  double next()
  {
    x += phi;
    x -= std::floor(x);
    return std::fma(2.0, x, -1.0);
  }
};

template <size_t n>
void subtest_solve_opt(utest::error_accumulator &ea)
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
    ea << utest::compare_numeric("wrong solve_opt<" + std::to_string(n) + ">", y[i], sum, eps<double> * std::sqrt(n));
  }
}

void test_solve_opt(utest::error_accumulator &ea)
{

  subtest_solve_opt<1>(ea);
  subtest_solve_opt<2>(ea);
  subtest_solve_opt<3>(ea);
  subtest_solve_opt<4>(ea);
  subtest_solve_opt<5>(ea);
  subtest_solve_opt<6>(ea);
  subtest_solve_opt<7>(ea);
  subtest_solve_opt<20>(ea);

  ea.throw_if_any();
}

void test_remove_tangent_components(utest::error_accumulator &ea)
{
  double u[2][3] = {
      {1.0, 0.0, 0.0},
      {0.0, std::sqrt(0.5), -std::sqrt(0.5)},
  };
  
  double x[3] = {0.741, 1.145, 1.876};
  
  const double rem = (x[1] + x[2]) / 2;

  math::remove_tangent_components(3, 2, x, u);
  
  ea << utest::compare_numeric("wrong remove_tangent_component 0", 0.0, x[0]);
  ea << utest::compare_numeric("wrong remove_tangent_component 1", rem, x[1]);
  ea << utest::compare_numeric("wrong remove_tangent_component 2", rem, x[2]);
}