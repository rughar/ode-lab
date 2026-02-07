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

  failed += utest::run(test_dot_product, "ling/dot_product");
  failed += utest::run(test_spectral_radius_estimate, "ling/spectral_radius_estimate");
  failed += utest::run(test_solve_opt, "ling/solve_opt");
  failed += utest::run(test_remove_tangent_components, "ling/remove_tangent_components");

  return failed;
}