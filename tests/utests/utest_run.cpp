#include <iostream>
#include <stdexcept>
#include <utest_frame.hpp>

#include <ling_test.hpp>

void test_smoke()
{
}

int main()
{
  int failed = 0;

  if (utest::run(test_smoke, "utest smoke"))
    return 1;

  failed += utest::run(test_dot_product, "dot_product");
  failed += utest::run(test_spectral_radius_estimate, "spectral_radius_estimate");

  return failed;
}