#include "driver/AsyncPipeline.hpp"

#include <thread>

namespace frehg2 {

void AsyncPipeline::advanceWindow(real t_begin, real dt) {
  const real t_end = t_begin + dt;

  if (!primed_) {
    // Priming step (k = 0): SW advances the first window, the exchange for it is applied, and
    // the groundwater catch-up is deferred so it can overlap the next surface window.
    cb_.swAdvance(t_end);
    cb_.applyExchange(dt);
    windows_.leading() = CouplingWindow{t_begin, t_end, dt, true};
    primed_ = true;
    return;
  }

  // Steady state (k >= 1): the trailing groundwater window is the previous leading window; its
  // exchange was already applied last call, so SW (new leading window) and GW (trailing window)
  // touch disjoint state and run concurrently. The exchange for the new leading window is
  // computed only after both join, exactly matching the sequential SW->exchange->GW ordering.
  const real gw_target = windows_.leading().t_end;
  std::thread gw_worker([this, gw_target]() { cb_.gwCatchUp(gw_target); });
  cb_.swAdvance(t_end);
  gw_worker.join();

  cb_.applyExchange(dt);
  windows_.swap();
  windows_.leading() = CouplingWindow{t_begin, t_end, dt, true};
}

void AsyncPipeline::drain() {
  if (!primed_) return;
  cb_.gwCatchUp(windows_.leading().t_end);
}

}  // namespace frehg2
