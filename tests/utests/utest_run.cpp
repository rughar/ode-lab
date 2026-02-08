#include <iostream>
#include <stdexcept>
#include <utest_frame.hpp>
#include <ling_test.hpp>

int main()
{
  utest::test_counter tc;

  utest::write_category("smoke");

  tc += utest::run(utest::test_smoke, "utest");
  if (tc.failed)
    return 1;

  utest::write_category("math::ling");

  tc += utest::run(test_dot_product, "dot_product");
  tc += utest::run(test_spectral_radius_estimate, "spectral_radius_estimate");
  tc += utest::run(test_solve_opt, "solve_opt");
  tc += utest::run(test_remove_tangent_components, "remove_tangent_components");

  return tc.failed ? 1 : 0;
}