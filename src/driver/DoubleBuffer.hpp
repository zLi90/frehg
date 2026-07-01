// Double buffer of coupling time windows for the async SW<->GW pipeline (P11.3.1).
//
// Architecture note (deviation from the P11 plan sketch, mirroring the P7 deviation
// recorded in Orchestrator.hpp): the plan sketch proposed two bundles each holding copies of
// the P2 `State`/`GwState` objects. In the realized P4/P5/P7 architecture each solver OWNS its
// own halo-padded fields (`SweFields` / `GwFields`), and the SWE and RE solvers are already two
// independent objects operating on disjoint data. The asynchronous overlap here is therefore
// NOT a snapshot of one domain's state into a second buffer — the surface solver never reads
// in-flight groundwater data (and vice versa); the only SW<->GW data motion is the exchange
// applied at the synchronization point. What actually has to be double-buffered is the
// *coupling time window* each domain is processing: the surface solver runs the LEADING window
// `[t_n, t_{n+1}]` while the groundwater solver catches up the TRAILING window `[t_{n-1}, t_n]`.
// This file models exactly that — two `CouplingWindow` slots with a `swap()` that rotates which
// slot is the leading (SW) window and which is the trailing (GW) window each step.
#ifndef FREHG2_DRIVER_DOUBLE_BUFFER_HPP
#define FREHG2_DRIVER_DOUBLE_BUFFER_HPP

#include "frehg2/core/define.hpp"

namespace frehg2 {

// One coupling window: the half-open simulated-time interval `[t_begin, t_end]` advanced by a
// single surface step, plus the surface step size `dt` used for the SW<->GW exchange volume.
struct CouplingWindow {
  real t_begin = 0.0;
  real t_end = 0.0;
  real dt = 0.0;
  bool valid = false;
};

// Two-slot double buffer. `leading()` is the window the surface solver is currently advancing;
// `trailing()` is the previous window the groundwater solver catches up concurrently. `swap()`
// rotates the two slots so the just-finished leading window becomes the trailing window for the
// next step (the classic ping-pong buffer, with no per-domain state copy required).
class DoubleBuffer {
 public:
  CouplingWindow& leading() { return flip_ ? b_ : a_; }
  CouplingWindow& trailing() { return flip_ ? a_ : b_; }
  const CouplingWindow& leading() const { return flip_ ? b_ : a_; }
  const CouplingWindow& trailing() const { return flip_ ? a_ : b_; }

  void swap() { flip_ = !flip_; }

 private:
  CouplingWindow a_;
  CouplingWindow b_;
  bool flip_ = false;
};

}  // namespace frehg2

#endif  // FREHG2_DRIVER_DOUBLE_BUFFER_HPP
