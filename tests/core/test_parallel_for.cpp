// P9 acceptance: the on-node parallel-loop helpers (parallelForRange / parallelForSurface /
// parallelForVolume) and Kokkos reductions produce results that are bit-identical to the plain
// serial loops they replace. This is the structural guarantee behind the Phase 9 conversion:
// every kernel was rewritten in terms of these wrappers, so if the wrappers reproduce serial
// element-wise writes and order-independent (min/max) reductions exactly, the kernels do too.
//
// These run on Kokkos::DefaultHostExecutionSpace (OpenMP on macOS, Serial fallback). The same
// code is GPU-capable; a real GPU execution check is deferred to the P10 Linux/NVIDIA validation
// (no CUDA on this machine), so the GPU equivalence case below is registered but skipped here.

#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include "frehg2_test.hpp"

#include <Kokkos_Core.hpp>

#include "core/Grid.hpp"
#include "frehg2/core/ParallelFor.hpp"
#include "frehg2/core/define.hpp"
#include "frehg2/core/types.hpp"

using namespace frehg2;

namespace {

// A reproducible, non-trivial per-cell value so that any reordering or skipped cell shows up as
// a mismatch rather than coincidentally agreeing.
KOKKOS_INLINE_FUNCTION real cellValue(int idx) {
  return 0.5 * static_cast<real>(idx) - 3.25 * static_cast<real>((idx % 7)) + 1.0;
}

}  // namespace

TEST_CASE("parallelForRange reproduces a serial whole-array write bit-for-bit") {
  const int n = 4096;
  RealArr1DHost serial("serial", n);
  RealArr1DHost par("par", n);
  for (int idx = 0; idx < n; ++idx) serial(idx) = cellValue(idx);
  parallelForRange("range", n, KOKKOS_LAMBDA(int idx) { par(idx) = cellValue(idx); });
  Kokkos::fence();
  for (int idx = 0; idx < n; ++idx) REQUIRE(par(idx) == serial(idx));
}

TEST_CASE("parallelForSurface reproduces a serial 2-D surface sweep bit-for-bit") {
  const int nx = 37, ny = 29;
  Grid g(nx, ny, 1, 1.0, 1.0, 1.0);
  RealArr1DHost serial("serial", g.nSurfaceStorageCell());
  RealArr1DHost par("par", g.nSurfaceStorageCell());
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) serial(g.getSurfaceIndex(i, j)) = cellValue(g.getSurfaceIndex(i, j));
  parallelForSurface("surf", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int idx = g.getSurfaceIndex(i, j);
    par(idx) = cellValue(idx);
  });
  Kokkos::fence();
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int idx = g.getSurfaceIndex(i, j);
      REQUIRE(par(idx) == serial(idx));
    }
}

TEST_CASE("parallelForVolume reproduces a serial 3-D volume sweep bit-for-bit") {
  const int nx = 11, ny = 13, nz = 7;
  Grid g(nx, ny, nz, 1.0, 1.0, 0.1, 1.0);
  RealArr1DHost serial("serial", g.nCell());
  RealArr1DHost par("par", g.nCell());
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) serial(g.getIndex(i, j, k)) = cellValue(g.getIndex(i, j, k));
  parallelForVolume("vol", nx, ny, nz, KOKKOS_LAMBDA(int i, int j, int k) {
    const int idx = g.getIndex(i, j, k);
    par(idx) = cellValue(idx);
  });
  Kokkos::fence();
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const int idx = g.getIndex(i, j, k);
        REQUIRE(par(idx) == serial(idx));
      }
}

TEST_CASE("parallel max-reduction is order-independent and matches the serial max") {
  const int n = 5000;
  RealArr1DHost a("a", n);
  real serial_max = -1.0e300;
  for (int idx = 0; idx < n; ++idx) {
    a(idx) = cellValue(idx);
    if (a(idx) > serial_max) serial_max = a(idx);
  }
  real par_max = -1.0e300;
  Kokkos::parallel_reduce(
      "max", Kokkos::RangePolicy<LoopExec>(0, n),
      KOKKOS_LAMBDA(int idx, real& m) {
        if (a(idx) > m) m = a(idx);
      },
      Kokkos::Max<real>(par_max));
  Kokkos::fence();
  // MAX is associative/commutative, so the floating-point result is exactly the serial scan.
  REQUIRE(par_max == serial_max);
}

TEST_CASE("parallel min-reduction is order-independent and matches the serial min") {
  const int n = 5000;
  RealArr1DHost a("a", n);
  real serial_min = 1.0e300;
  for (int idx = 0; idx < n; ++idx) {
    a(idx) = cellValue(idx);
    if (a(idx) < serial_min) serial_min = a(idx);
  }
  real par_min = 1.0e300;
  Kokkos::parallel_reduce(
      "min", Kokkos::RangePolicy<LoopExec>(0, n),
      KOKKOS_LAMBDA(int idx, real& m) {
        if (a(idx) < m) m = a(idx);
      },
      Kokkos::Min<real>(par_min));
  Kokkos::fence();
  REQUIRE(par_min == serial_min);
}

TEST_CASE("parallel skip-predicate (early return) matches a serial continue") {
  // The P9 conversions translate `if (!active) continue;` into `if (!active) return;` inside the
  // lambda. Verify the two skip semantics leave identical results (untouched cells keep their
  // initial value; touched cells get the computed value).
  const int n = 2048;
  RealArr1DHost serial("serial", n);
  RealArr1DHost par("par", n);
  Kokkos::deep_copy(serial, -1.0);
  Kokkos::deep_copy(par, -1.0);
  auto activeCell = [](int idx) { return (idx % 3) != 0; };
  for (int idx = 0; idx < n; ++idx) {
    if (!activeCell(idx)) continue;
    serial(idx) = cellValue(idx);
  }
  parallelForRange("skip", n, KOKKOS_LAMBDA(int idx) {
    if ((idx % 3) == 0) return;
    par(idx) = cellValue(idx);
  });
  Kokkos::fence();
  for (int idx = 0; idx < n; ++idx) REQUIRE(par(idx) == serial(idx));
}

// GPU host/device equivalence is registered as a separate CTest entry (label "gpu", DISABLED on
// this machine) in tests/core/CMakeLists.txt — there is no CUDA here, so it is skipped per
// docs/gpu_validation_policy.md and deferred to the P10 Linux/NVIDIA validation. On a Kokkos-CUDA
// build the same wrappers (LoopExec → Cuda) run on the device and are compared to the host
// reference above.
