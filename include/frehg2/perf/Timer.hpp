// Kokkos-backed scoped timer with optional Kokkos profiling regions (P21 Task 21.3.1).
#ifndef FREHG2_PERF_TIMER_HPP
#define FREHG2_PERF_TIMER_HPP

#include <Kokkos_Core.hpp>

#include "frehg2/perf/PerfRegions.hpp"

namespace frehg2 {

class PerfRecorder;

namespace perf {

// RAII timer: accumulates elapsed seconds into `recorder` on destruction and emits a
// Kokkos::Profiling region for nvprof / vtune / Instruments when profiling is enabled.
class ScopedTimer {
 public:
  ScopedTimer(PerfRecorder* recorder, Region region);
  ~ScopedTimer();

  ScopedTimer(const ScopedTimer&) = delete;
  ScopedTimer& operator=(const ScopedTimer&) = delete;

 private:
  PerfRecorder* recorder_;
  Region region_;
  Kokkos::Timer timer_;
  bool active_;
};

}  // namespace perf
}  // namespace frehg2

#endif  // FREHG2_PERF_TIMER_HPP
