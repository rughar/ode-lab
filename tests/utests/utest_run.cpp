#include <iostream>
#include <stdexcept>
#include <utest_frame.hpp>

void test_dummy()
{
  throw std::runtime_error("Newton solver did not converge");
}

int main() 
{
  int failed = 0;
  failed += utest::run(test_dummy, "Dummy test"); 
  return failed;
}