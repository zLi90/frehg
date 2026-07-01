// Asynchronous SW<->GW coupling pipeline (P11.3.2) — CPU/OpenMP realization.
//
// Promotes the P6 blocking-sync coupling to a double-buffered pipeline: while the surface
// solver advances the LEADING window `[t_n, t_{n+1}]`, the groundwater solver catches up the
// TRAILING window `[t_{n-1}, t_n]` on a worker thread. The coupling flux is computed at the
// synchronization point (after both threads join) from the converged SW and GW states, exactly
// as in the synchronous path — so the async path executes the identical Gauss–Seidel operation
// sequence (SW advance -> exchange -> GW catch-up) and is numerically bit-identical to the
// sequential coupling, not merely within 1e-10.
//
// Threading / PETSc note: the two solver advances dispatch to two OS threads (this is the
// "thread pool with per-domain locks" the plan calls for), but the actual `KSPSolve` calls are
// serialized by a caller-supplied mutex because PETSc `KSP` solves on a shared communicator are
// not thread-safe. On a single rank (the validated CPU path) the surface and groundwater solves
// therefore alternate under the lock while their independent host work overlaps; true
// concurrent solves require the per-domain communicator split (`PetscSubcomm` + CUDA/HIP
// streams) deferred to the GPU pass (P10, Task 10.3.7). No GPU-bridge symbol is referenced
// here: P11 is CPU/OpenMP only.
//
// No PETSc type appears in this header or its translation unit: the pipeline drives the domains
// only through the three callbacks the Orchestrator binds (advance SW, catch up GW, apply
// exchange), each of which already routes through the backend-agnostic LinearSolver seam.
#ifndef FREHG2_DRIVER_ASYNC_PIPELINE_HPP
#define FREHG2_DRIVER_ASYNC_PIPELINE_HPP

#include <functional>

#include "driver/DoubleBuffer.hpp"
#include "frehg2/core/define.hpp"

namespace frehg2 {

class AsyncPipeline {
 public:
  // Callbacks bound by the Orchestrator to its private solver operations.
  struct Callbacks {
    // Advance the surface solver by ONE surface step so it reaches simulated time `t_target`
    // (applies rainfall/evaporation for that step). Must take the PETSc-solve lock internally.
    std::function<void(real t_target)> swAdvance;
    // Advance the groundwater solver (one or more adaptive substeps) until it reaches
    // `t_target`. Must take the PETSc-solve lock around each substep solve internally.
    std::function<void(real t_target)> gwCatchUp;
    // Compute the SW<->GW exchange flux from the converged SW and GW states and apply it to
    // both domains for surface step size `dt` (no solver advance). Runs on the calling thread
    // after the join, so it needs no lock.
    std::function<void(real dt)> applyExchange;
  };

  explicit AsyncPipeline(Callbacks cb) : cb_(std::move(cb)) {}

  // Process the surface window `[t_begin, t_begin + dt]`.
  //
  // First (priming) call: advance SW over the window and apply its exchange; the groundwater
  // catch-up for this window is deferred (it overlaps the NEXT surface window). Subsequent
  // calls: launch the groundwater catch-up of the previous (trailing) window on a worker
  // thread, advance the surface solver over the new leading window on the calling thread, join,
  // then compute+apply the exchange for the new leading window. After this returns the surface
  // solver is at `t_begin + dt` and the groundwater solver is one window behind (at `t_begin`).
  void advanceWindow(real t_begin, real dt);

  // Finish the pending trailing groundwater window so the groundwater solver is synchronized
  // with the surface solver (both at the leading window's end). Idempotent: a second call with
  // no intervening advanceWindow() is a no-op because the catch-up target is already reached.
  // Must be called before any output/checkpoint that reads both domains, and once at run end.
  void drain();

  bool primed() const { return primed_; }

 private:
  Callbacks cb_;
  DoubleBuffer windows_;
  bool primed_ = false;
};

}  // namespace frehg2

#endif  // FREHG2_DRIVER_ASYNC_PIPELINE_HPP
