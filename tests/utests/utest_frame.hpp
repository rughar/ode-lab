#pragma once

#include <stdexcept>
#include <string>
#include <iostream>
#include <format>
#include <string_view>

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

    ~test_counter()
    {
      std::cout << "\033[33mTest summary\033[0m: " << total - failed << " passed, " << failed << " failed, out of " << total << " tests.\n";
    }
  };

  class error_accumulator
  {
  private:
    std::string msg;
    std::string p_prefix = "\033[32mPASSED\033[0m: ";
    std::string f_prefix = "\033[31mFAILED\033[0m: ";
    bool silent;

  public:
    error_accumulator(const bool silence) : silent(silence)
    {
    }

    error_accumulator &operator<<(const std::string_view &error_message)
    {
      if (!error_message.empty())
        msg += "\n  - " + std::string(error_message);
      return *this;
    }

    bool is_silent() const
    {
      return silent;
    }

    void throw_if_any() const
    {
      if (!msg.empty())
        throw std::runtime_error(msg);
    }
  };

  template <class TestFn>
  int run(TestFn &&test_func, const std::string_view &test_name) noexcept
  {
    error_accumulator ea(test_name.empty());

    const std::string passed_prefix = "\033[32mPASSED\033[0m: " + std::string(test_name);
    const std::string failed_prefix = "\033[31mFAILED\033[0m: " + std::string(test_name);

    try
    {
      test_func(ea);
      ea.throw_if_any();
      if (!ea.is_silent())
        std::cout << passed_prefix << "\n";
      return 0;
    }
    catch (const std::exception &e)
    {
      if (!ea.is_silent())
        std::cout << failed_prefix << e.what() << "\n";
      return 1;
    }
    catch (...)
    {
      if (!ea.is_silent())
        std::cout << failed_prefix << " (Unknown exception)\n";
      return 1;
    }
  }

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

  template <class U>
  std::string compare_numeric(const std::string_view &msg, U expected, U actual, U tol = U(0))
  {
    if (std::abs(expected - actual) > tol)
    {
      std::string expected_str = std::format("{:.17g}", expected);
      std::string actual_str = std::format("{:.17g}", actual);
      highlight_difference(expected_str, actual_str);
      return std::string(msg) + ": expected " + expected_str + ", got " + actual_str;
    }
    return "";
  }

  void write_category(const std::string_view &category_name)
  {
    std::cout << "\033[34m" << category_name << "\033[0m\n";
  }

  void test_smoke(error_accumulator &ea)
  {
    if (1 + 1 != 2)
      ea << "math is broken";

    auto must_fail = [&ea](error_accumulator &ea_in)
    {
      ea_in << "This test should failed silently. Becouse you are seeing this message, it did not.";
      if (!ea_in.is_silent())
        ea << "Silencing did not work in the subtest.";
    };

    if (utest::run(must_fail, "") == 0)
      ea << "run() did not catch the exception";
  }
};