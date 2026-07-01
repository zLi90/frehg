#include "ic/FormulaParser.hpp"

#include <cctype>
#include <cmath>
#include <stdexcept>
#include <string>

namespace frehg2 {

namespace {

bool isIdentStart(char c) { return std::isalpha(static_cast<unsigned char>(c)) || c == '_'; }
bool isIdentChar(char c) {
  return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
}

}  // namespace

int FormulaParser::opPrecedence(TokKind k) {
  switch (k) {
    case TokKind::Caret:
      return 4;
    case TokKind::Star:
    case TokKind::Slash:
      return 3;
    case TokKind::Plus:
    case TokKind::Minus:
      return 2;
    default:
      return 0;
  }
}

bool FormulaParser::isRightAssociative(TokKind k) { return k == TokKind::Caret; }

bool FormulaParser::isFuncToken(TokKind k) {
  return k == TokKind::Sin || k == TokKind::Cos || k == TokKind::Exp || k == TokKind::Sqrt ||
         k == TokKind::Abs;
}

FormulaParser::FormulaParser(std::string expr) : expr_(std::move(expr)) { compile(); }

void FormulaParser::compile() {
  std::vector<Token> output;
  std::vector<Token> ops;
  const std::string& s = expr_;
  size_t i = 0;
  bool expect_operand = true;

  auto push_op = [&](Token t) {
    while (!ops.empty()) {
      const Token top = ops.back();
      if (top.kind == TokKind::LParen) break;
      if (opPrecedence(top.kind) > opPrecedence(t.kind) ||
          (opPrecedence(top.kind) == opPrecedence(t.kind) && !isRightAssociative(t.kind))) {
        output.push_back(top);
        ops.pop_back();
      } else {
        break;
      }
    }
    ops.push_back(t);
  };

  while (i < s.size()) {
    const char c = s[i];
    if (std::isspace(static_cast<unsigned char>(c))) {
      ++i;
      continue;
    }
    if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
      size_t j = i + 1;
      while (j < s.size()) {
        const char d = s[j];
        if (std::isdigit(static_cast<unsigned char>(d)) || d == '.') {
          ++j;
          continue;
        }
        if (d == 'e' || d == 'E') {
          ++j;
          if (j < s.size() && (s[j] == '+' || s[j] == '-')) ++j;
          continue;
        }
        break;
      }
      try {
        const real v = std::stod(s.substr(i, j - i));
        output.push_back({TokKind::Number, v});
      } catch (const std::exception&) {
        throw std::runtime_error("FormulaParser: invalid number near '" + s.substr(i, j - i) +
                                 "' in '" + expr_ + "'");
      }
      i = j;
      expect_operand = false;
      continue;
    }
    if (isIdentStart(c)) {
      size_t j = i + 1;
      while (j < s.size() && isIdentChar(s[j])) ++j;
      const std::string id = s.substr(i, j - i);
      if (id == "x")
        output.push_back({TokKind::VarX});
      else if (id == "y")
        output.push_back({TokKind::VarY});
      else if (id == "z")
        output.push_back({TokKind::VarZ});
      else if (id == "t")
        output.push_back({TokKind::VarT});
      else if (id == "pi")
        output.push_back({TokKind::Pi});
      else if (id == "e")
        output.push_back({TokKind::E});
      else if (id == "sin")
        ops.push_back({TokKind::Sin});
      else if (id == "cos")
        ops.push_back({TokKind::Cos});
      else if (id == "exp")
        ops.push_back({TokKind::Exp});
      else if (id == "sqrt")
        ops.push_back({TokKind::Sqrt});
      else if (id == "abs")
        ops.push_back({TokKind::Abs});
      else
        throw std::runtime_error("FormulaParser: unknown identifier '" + id + "' in '" + expr_ +
                                 "'");
      i = j;
      expect_operand = false;
      continue;
    }
    if (c == '(') {
      ops.push_back({TokKind::LParen});
      ++i;
      expect_operand = true;
      continue;
    }
    if (c == ')') {
      while (!ops.empty() && ops.back().kind != TokKind::LParen) {
        output.push_back(ops.back());
        ops.pop_back();
      }
      if (ops.empty())
        throw std::runtime_error("FormulaParser: unmatched ')' in '" + expr_ + "'");
      ops.pop_back();
      if (!ops.empty() && isFuncToken(ops.back().kind)) {
        output.push_back(ops.back());
        ops.pop_back();
      }
      ++i;
      expect_operand = false;
      continue;
    }
    TokKind op = TokKind::Plus;
    if (c == '+')
      op = TokKind::Plus;
    else if (c == '-')
      op = TokKind::Minus;
    else if (c == '*')
      op = TokKind::Star;
    else if (c == '/')
      op = TokKind::Slash;
    else if (c == '^')
      op = TokKind::Caret;
    else
      throw std::runtime_error("FormulaParser: unexpected character '" + std::string(1, c) +
                               "' in '" + expr_ + "'");

    if (expect_operand && op == TokKind::Minus) {
      ops.push_back({TokKind::Neg});
    } else if (expect_operand && op == TokKind::Plus) {
      // unary plus is a no-op
    } else {
      push_op({op});
      expect_operand = true;
    }
    ++i;
  }
  while (!ops.empty()) {
    if (ops.back().kind == TokKind::LParen)
      throw std::runtime_error("FormulaParser: unmatched '(' in '" + expr_ + "'");
    output.push_back(ops.back());
    ops.pop_back();
  }
  rpn_ = std::move(output);
}

real FormulaParser::eval(const FormulaVars& v) const {
  std::vector<real> st;
  st.reserve(rpn_.size());
  auto pop = [&]() -> real {
    if (st.empty()) throw std::runtime_error("FormulaParser: stack underflow in '" + expr_ + "'");
    const real x = st.back();
    st.pop_back();
    return x;
  };
  for (const Token& t : rpn_) {
    switch (t.kind) {
      case TokKind::Number:
        st.push_back(t.number);
        break;
      case TokKind::VarX:
        st.push_back(v.x);
        break;
      case TokKind::VarY:
        st.push_back(v.y);
        break;
      case TokKind::VarZ:
        st.push_back(v.z);
        break;
      case TokKind::VarT:
        st.push_back(v.t);
        break;
      case TokKind::Pi:
        st.push_back(static_cast<real>(M_PI));
        break;
      case TokKind::E:
        st.push_back(static_cast<real>(M_E));
        break;
      case TokKind::Neg: {
        const real a = pop();
        st.push_back(-a);
        break;
      }
      case TokKind::Plus: {
        const real b = pop(), a = pop();
        st.push_back(a + b);
        break;
      }
      case TokKind::Minus: {
        const real b = pop(), a = pop();
        st.push_back(a - b);
        break;
      }
      case TokKind::Star: {
        const real b = pop(), a = pop();
        st.push_back(a * b);
        break;
      }
      case TokKind::Slash: {
        const real b = pop(), a = pop();
        st.push_back(a / b);
        break;
      }
      case TokKind::Caret: {
        const real b = pop(), a = pop();
        st.push_back(std::pow(a, b));
        break;
      }
      case TokKind::Sin: {
        const real a = pop();
        st.push_back(std::sin(a));
        break;
      }
      case TokKind::Cos: {
        const real a = pop();
        st.push_back(std::cos(a));
        break;
      }
      case TokKind::Exp: {
        const real a = pop();
        st.push_back(std::exp(a));
        break;
      }
      case TokKind::Sqrt: {
        const real a = pop();
        st.push_back(std::sqrt(a));
        break;
      }
      case TokKind::Abs: {
        const real a = pop();
        st.push_back(std::fabs(a));
        break;
      }
      default:
        break;
    }
  }
  if (st.size() != 1)
    throw std::runtime_error("FormulaParser: invalid expression '" + expr_ + "' (stack=" +
                             std::to_string(st.size()) + ", rpn=" + std::to_string(rpn_.size()) +
                             ")");
  return st[0];
}

}  // namespace frehg2
