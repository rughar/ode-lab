#include <iostream>
#include <utest_frame.hpp>

void test_dummy()
{
}

int main() 
{
  int failed = 0;
  failed += utest::run(test_dummy, "Dummy test"); 
  return failed;
}