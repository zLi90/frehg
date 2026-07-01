#include "frehg2/perf/Timer.hpp"

#include "frehg2/perf/PerfRecorder.hpp"

namespace frehg2 {
namespace perf {

ScopedTimer::ScopedTimer(PerfRecorder* recorder, Region region)
    : recorder_(recorder), region_(region), active_(recorder != nullptr) {
  if (!active_) return;
#if defined(KOKKOS_ENABLE_PROFILING)
  Kokkos::Profiling::pushRegion(regionName(region_));
#else
  (void)region_;
#endif
}

ScopedTimer::~ScopedTimer() {
  if (!active_) return;
  recorder_->addTime(region_, timer_.seconds());
#if defined(KOKKOS_ENABLE_PROFILING)
  Kokkos::Profiling::popRegion();
#endif
}

}  // namespace perf
}  // namespace frehg2
