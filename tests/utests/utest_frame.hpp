#pragma once

#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <format>

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

  struct error_accumulator
  {
    std::vector<std::string> errors;

    void add(const std::string &error)
    {
      errors.push_back(error);
    }

    void throw_if_any()
    {
      if (!errors.empty())
      {
        std::string msg;
        std::string sep;

        for (const auto &error : errors) {
          msg += sep + error;
          sep = "\n  - ";
        }
        throw std::runtime_error(msg);
      }
    }
  };

  void highlight_difference(std::string &a, std::string &b)
  {
    const std::size_t n = std::min(a.size(), b.size());

    std::size_t diff = n;
    for (std::size_t i = 0; i < n; ++i)
    {
      if (a[i] != b[i])
      {
        diff = i;
        break;
      }
    }

    if (a.size() != b.size())
      diff = std::min(diff, n);

    if (diff == n && a.size() == b.size())
      return;

    a.insert(diff, "\033[38;5;208m");
    a.append("\033[0m");

    b.insert(diff, "\033[38;5;208m");
    b.append("\033[0m");
  }

  template <class S, class U>
  void compare_numeric(const S &msg, U expected, U actual, U tol = std::numeric_limits<U>::epsilon())
  {
    if (std::abs(expected - actual) > tol)
    {
      std::string expected_str = std::format("{:.17g}", expected);
      std::string actual_str = std::format("{:.17g}", actual);
      highlight_difference(expected_str, actual_str);
      throw std::runtime_error(std::string(msg) + ": expected " + expected_str + ", got " + actual_str);
    }
  }

  template <class S>
  void write_category(const S &category_name)
  {
    std::cout << "\033[34m" << category_name << "\033[0m\n";
  }

  template <class S, class TestFn>
  int run(TestFn &&test_func, const S &test_name) noexcept
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
        std::cout << "\033[31mFAILED\033[0m: " << test_name << "\n  - " << e.what() << "\n";
      return 1;
    }
    catch (...)
    {
      if (*test_name)
        std::cout << "\033[31mFAILED\033[0m: " << test_name << "\n  - " << "Unknown exception\n";
      return 1;
    }
  }

  void test_smoke()
  {
    if (1 + 1 != 2)
      throw std::runtime_error("math is broken");

    auto must_throw = []()
    { throw std::runtime_error("expected failure"); };

    if (utest::run(must_throw, "") == 0)
      throw std::runtime_error("run() did not catch the exception");
  }
};