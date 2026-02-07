#pragma once

#include <stdexcept>

namespace utest
{
  struct test_counter
  {
    size_t total = 0;
    size_t failed = 0;

    void operator+=(size_t result)
    {
      ++total;
      if (result != 0)
        ++failed;
    };

    void write_summary() const
    {
      std::cout << "\033[33mTest summary\033[0m: " << total - failed << " passed, " << failed << " failed, out of " << total << " tests.\n";
    }
  };

  void write_category(const char *category_name)
  {
    std::cout << "\033[34m" << category_name << "\033[0m\n";
  }

  template <class TestFn>
  int run(TestFn &&test_func, const char *test_name) noexcept
  {
    try
    {
      test_func();
      if (*test_name)
        std::cout << "\033[32mPASSED\033[0m: " << test_name << "\n";
      return 0;
    }
    catch (const std::exception &e)
    {
      if (*test_name)
        std::cerr << "\033[31mFAILED\033[0m: " << test_name << " - " << e.what() << "\n";
      return 1;
    }
    catch (...)
    {
      if (*test_name)
        std::cerr << "\033[31mFAILED\033[0m: " << test_name << " - " << "Unknown exception\n";
      return 1;
    }
  }
};