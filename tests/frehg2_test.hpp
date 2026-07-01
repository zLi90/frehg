// Minimal Catch2-style test harness (vendored) for Frehg2.
//
// Rationale (documented per .cursorrules "Catch2 v3 preferred, ... vendored ... allowed
// if P1 documents the choice"): no Catch2/GTest is installed in the local prefix and the
// build must stay network-free and deterministic. This header provides a tiny subset of
// the Catch2 API (TEST_CASE / REQUIRE / CHECK / REQUIRE_FALSE / Approx) so test sources
// are Catch2-compatible and can migrate to real Catch2 v3 later without edits.
//
// Usage (one executable per test file):
//   #define FREHG2_TEST_IMPL          // exactly one TU provides main()
//   #define FREHG2_TEST_USE_KOKKOS    // optional: wrap run in Kokkos init/finalize
//   #include "frehg2_test.hpp"
//   TEST_CASE("name") { REQUIRE(x == Approx(y).margin(1e-9)); }
#ifndef FREHG2_TEST_HPP
#define FREHG2_TEST_HPP

#include <cmath>
#include <cstdio>
#include <exception>
#include <functional>
#include <string>
#include <vector>

namespace frehg2test {

struct TestCase {
  std::string name;
  std::function<void()> fn;
};

inline std::vector<TestCase>& registry() {
  static std::vector<TestCase> r;
  return r;
}

inline int& assertionFailures() {
  static int n = 0;
  return n;
}

// Thrown by REQUIRE on failure to abort the current test case.
struct AssertionFailure : public std::exception {
  const char* what() const noexcept override { return "REQUIRE failed"; }
};

struct Registrar {
  Registrar(const std::string& name, std::function<void()> fn) {
    registry().push_back({name, std::move(fn)});
  }
};

// Catch2-style approximate comparison: passes if within an absolute margin OR a relative
// epsilon of the reference value.
class Approx {
 public:
  explicit Approx(double value) : value_(value) {}
  Approx& margin(double m) {
    margin_ = m;
    return *this;
  }
  Approx& epsilon(double e) {
    epsilon_ = e;
    return *this;
  }
  bool matches(double other) const {
    const double diff = std::fabs(other - value_);
    if (diff <= margin_) {
      return true;
    }
    const double scale = std::fmax(std::fabs(value_), std::fabs(other));
    return diff <= epsilon_ * scale;
  }

 private:
  double value_;
  double margin_ = 0.0;
  double epsilon_ = 1e-12;
};

inline bool operator==(double lhs, const Approx& rhs) { return rhs.matches(lhs); }
inline bool operator==(const Approx& lhs, double rhs) { return lhs.matches(rhs); }
inline bool operator!=(double lhs, const Approx& rhs) { return !rhs.matches(lhs); }
inline bool operator!=(const Approx& lhs, double rhs) { return !lhs.matches(rhs); }

inline void reportFailure(const char* expr, const char* file, int line, bool fatal) {
  std::fprintf(stderr, "  %s assertion failed: %s\n    at %s:%d\n",
               fatal ? "REQUIRE" : "CHECK", expr, file, line);
  ++assertionFailures();
}

inline int runAll() {
  int failed_cases = 0;
  for (const auto& tc : registry()) {
    const int before = assertionFailures();
    bool aborted = false;
    try {
      tc.fn();
    } catch (const AssertionFailure&) {
      aborted = true;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "  unexpected exception in '%s': %s\n", tc.name.c_str(),
                   e.what());
      ++assertionFailures();
    }
    const bool ok = (assertionFailures() == before) && !aborted;
    std::fprintf(stderr, "[%s] %s\n", ok ? "PASS" : "FAIL", tc.name.c_str());
    if (!ok) {
      ++failed_cases;
    }
  }
  std::fprintf(stderr, "\n%d test case(s), %d failed, %d assertion failure(s)\n",
               static_cast<int>(registry().size()), failed_cases, assertionFailures());
  return failed_cases == 0 ? 0 : 1;
}

}  // namespace frehg2test

#define FREHG2_CONCAT_INNER(a, b) a##b
#define FREHG2_CONCAT(a, b) FREHG2_CONCAT_INNER(a, b)

#define TEST_CASE(name)                                                       \
  static void FREHG2_CONCAT(frehg2_test_fn_, __LINE__)();                     \
  static ::frehg2test::Registrar FREHG2_CONCAT(frehg2_test_reg_, __LINE__)(   \
      name, &FREHG2_CONCAT(frehg2_test_fn_, __LINE__));                       \
  static void FREHG2_CONCAT(frehg2_test_fn_, __LINE__)()

#define REQUIRE(expr)                                                  \
  do {                                                                 \
    if (!(expr)) {                                                     \
      ::frehg2test::reportFailure(#expr, __FILE__, __LINE__, true);    \
      throw ::frehg2test::AssertionFailure();                         \
    }                                                                  \
  } while (false)

#define REQUIRE_FALSE(expr) REQUIRE(!(expr))

#define CHECK(expr)                                                    \
  do {                                                                 \
    if (!(expr)) {                                                     \
      ::frehg2test::reportFailure(#expr, __FILE__, __LINE__, false);   \
    }                                                                  \
  } while (false)

using ::frehg2test::Approx;

#ifdef FREHG2_TEST_IMPL
#ifdef FREHG2_TEST_USE_KOKKOS
#include <Kokkos_Core.hpp>
int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  const int rc = ::frehg2test::runAll();
  Kokkos::finalize();
  return rc;
}
#else
int main(int argc, char** argv) {
  (void)argc;
  (void)argv;
  return ::frehg2test::runAll();
}
#endif  // FREHG2_TEST_USE_KOKKOS
#endif  // FREHG2_TEST_IMPL

#endif  // FREHG2_TEST_HPP
