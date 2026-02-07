#include <iostream>
#include <stdexcept>
#include <vector>
#include <utest_frame.hpp>
#include <ling_test.hpp>

void test_smoke()
{
  if (1 + 1 != 2)
    throw std::runtime_error("math is broken");

  auto must_throw = []()
  { throw std::runtime_error("expected failure"); };
  
  if (utest::run(must_throw, "") == 0)
    throw std::runtime_error("run() did not catch the exception");
}

int main()
{
  utest::test_counter tc;
  utest::write_category("smoke");

  tc += utest::run(test_smoke, "utest smoke");
  if (tc.failed)
    return 1;

  utest::write_category("math::ling");

  tc += utest::run(test_dot_product, "dot_product");
  tc += utest::run(test_spectral_radius_estimate, "spectral_radius_estimate");
  tc += utest::run(test_solve_opt, "solve_opt");
  tc += utest::run(test_remove_tangent_components, "remove_tangent_components");

  tc.write_summary();

  return tc.failed ? 1 : 0;
}