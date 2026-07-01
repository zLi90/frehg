// CQ1 negative fixture (P1.5). This file is DELIBERATELY not warning-clean: it contains
// an unused local variable, which -Wall promotes to an error under -Werror.
//
// It is NOT part of any add_executable target. The CTest test
// `cq1_strict_warnings_rejects_bad_code` compiles this file with the strict warning flags
// and is marked WILL_FAIL, so the suite passes only if strict mode correctly REJECTS it.
// Note: -Wno-unused-parameter is in the strict set, so the trigger is an unused *local*,
// not an unused parameter.
int main() {
  int unused_local = 42;  // -Wunused-variable / -Wunused-but-set-variable -> error
  return 0;
}
