// Per-run performance counters (P21 Task 21.3.2). Per-rank only — no global aggregation.
#ifndef FREHG2_PERF_COUNTERS_HPP
#define FREHG2_PERF_COUNTERS_HPP

#include <cstdint>

namespace frehg2 {
namespace perf {

struct Counters {
  int64_t cells_touched = 0;   // owned cells advanced this run (SW + GW accumulators)
  int64_t ksp_iterations = 0;  // sum of LinearSolver::getIterationCount() per solve
  int64_t bytes_staged = 0;    // approximate COO/RHS bytes moved host↔solver per assembly

  void addCells(int64_t n) { cells_touched += n; }
  void addKspIters(int n) { ksp_iterations += n; }
  void addBytes(int64_t n) { bytes_staged += n; }
};

}  // namespace perf
}  // namespace frehg2

#endif  // FREHG2_PERF_COUNTERS_HPP
