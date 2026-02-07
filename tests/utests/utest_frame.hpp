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
      std::cout << "\033[32mPASSED\033[0m: " << test_name << "\n";
      return 0;
    }
    catch (const std::exception &e)
    {
      std::cerr << "\033[31mFAILED\033[0m: " << test_name << " - " << e.what() << "\n";
      return 1;
    }
    catch (...)
    {
      std::cerr << "\033[31mFAILED\033[0m: " << test_name << " - " << "Unknown exception\n";
      return 1;
    }
  }
};