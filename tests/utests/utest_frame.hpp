#pragma once

#include <stdexcept>

namespace utest
{
  template <class TestFn>
  int run(TestFn &&test_func, const char *test_name) noexcept
  {
    try
    {
      test_func();
      std::cout << "Passed: " << test_name << "\n";
      return 0;
    }
    catch (const std::exception &e)
    {
      std::cerr << "Failed: " << test_name << " - " << e.what() << "\n";
      return 1;
    }
    catch (...)
    {
      std::cerr << "Failed: " << test_name << " - " << "Unknown exception\n";
      return 1;
    }
  }
};