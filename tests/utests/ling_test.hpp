#pragma once
#include <utest_frame.hpp>
#include <ling.hpp>

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
    1.0, 2.0,
    0.0, 1.0,
  };
  double expected = 1.0;
  double result = math::spectral_radius_estimate(2, mat);
  if (result != expected)
    throw std::runtime_error("wrong spectral radius on 2*2 matrix");
}