// P14.3.2: formula parser correctness on 10+ expressions (shunting-yard RPN).
#define FREHG2_TEST_IMPL
#include "frehg2_test.hpp"

#include <cmath>

#include "ic/FormulaParser.hpp"

using namespace frehg2;

namespace {

constexpr double kTol = 1.0e-12;

void checkExpr(const char* expr, const FormulaVars& v, double expected) {
  FormulaParser p(expr);
  const double got = p.eval(v);
  REQUIRE(got == Approx(expected).margin(kTol));
}

}  // namespace

TEST_CASE("formula parser literals and arithmetic") {
  FormulaVars v;
  checkExpr("0", v, 0.0);
  checkExpr("0.5+0.5", v, 1.0);
  checkExpr("2*3", v, 6.0);
  checkExpr("1+2*3", v, 7.0);
  checkExpr("(1+2)*3", v, 9.0);
  checkExpr("2^3", v, 8.0);
  checkExpr("10-3*2", v, 4.0);
  checkExpr("10/4", v, 2.5);
}

TEST_CASE("formula parser constants and unary minus") {
  FormulaVars v;
  checkExpr("pi", v, M_PI);
  checkExpr("e", v, M_E);
  v.x = 2.0;
  checkExpr("-x", v, -2.0);
  checkExpr("+x", v, 2.0);
}

TEST_CASE("formula parser functions") {
  FormulaVars v;
  checkExpr("sin(0)", v, 0.0);
  checkExpr("cos(0)", v, 1.0);
  checkExpr("sqrt(4)", v, 2.0);
  checkExpr("exp(0)", v, 1.0);
  checkExpr("abs(-3.5)", v, 3.5);
  v.x = 0.5;
  checkExpr("sin(pi*x)", v, 1.0);
}

TEST_CASE("formula parser spatial Gaussian") {
  FormulaVars v;
  v.x = 0.0;
  v.y = 0.0;
  checkExpr("exp(-(x^2+y^2)/0.01)", v, 1.0);
  v.x = 1.0;
  v.y = 0.0;
  checkExpr("exp(-(x^2+y^2)/0.01)", v, std::exp(-100.0));
}

TEST_CASE("formula parser uses z and t") {
  FormulaVars v;
  v.z = 2.0;
  v.t = 3.0;
  checkExpr("z+t", v, 5.0);
}
