// Tiny expression parser for spatial IC formulas (P14.3.2).
//
// Tokenizes and evaluates expressions limited to `+ - * / ^ ( )`, numeric literals, identifiers
// `x y z t pi e`, and unary/binary functions `sin cos exp sqrt abs`. Uses shunting-yard → RPN
// on the host (ICs are applied once at init, not inside Kokkos lambdas). Throws on unknown
// tokens or mismatched parentheses.
#ifndef FREHG2_IC_FORMULA_PARSER_HPP
#define FREHG2_IC_FORMULA_PARSER_HPP

#include <string>
#include <vector>

#include "frehg2/core/define.hpp"

namespace frehg2 {

struct FormulaVars {
  real x = 0.0;
  real y = 0.0;
  real z = 0.0;
  real t = 0.0;
};

class FormulaParser {
 public:
  explicit FormulaParser(std::string expr);

  // Evaluate the compiled RPN at the given variable values.
  real eval(const FormulaVars& v) const;

  const std::string& expression() const { return expr_; }

 private:
  enum class TokKind {
    Number,
    VarX,
    VarY,
    VarZ,
    VarT,
    Pi,
    E,
    Plus,
    Minus,
    Star,
    Slash,
    Caret,
    LParen,
    RParen,
    Sin,
    Cos,
    Exp,
    Sqrt,
    Abs,
    Neg
  };

  struct Token {
    TokKind kind;
    real number = 0.0;
  };

  std::string expr_;
  std::vector<Token> rpn_;

  static int opPrecedence(TokKind k);
  static bool isRightAssociative(TokKind k);
  static bool isFuncToken(TokKind k);

  void compile();
};

}  // namespace frehg2

#endif  // FREHG2_IC_FORMULA_PARSER_HPP
