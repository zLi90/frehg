// P2.2 acceptance: State / GwState dimensions, deep-copy isolation, time-level swap.
#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include "frehg2_test.hpp"

#include "core/GwState.hpp"
#include "core/State.hpp"

using namespace frehg2;

TEST_CASE("surface state dimensions match grid") {
  Grid g(8, 5, 1, 1.0, 1.0, 1.0);
  State s(g);
  REQUIRE(s.eta.extent(0) == 40u);
  REQUIRE(s.u.extent(0) == 40u);
  REQUIRE(s.v.extent(0) == 40u);
  REQUIRE(s.h.extent(0) == 40u);
  REQUIRE(s.qss.extent(0) == 40u);
}

TEST_CASE("gw state dimensions match active cells") {
  Grid g(4, 3, 5, 1.0, 1.0, 0.1, 1.0);
  GwState s(g);
  REQUIRE(s.h.extent(0) == 60u);
  REQUIRE(s.wc.extent(0) == 60u);
  REQUIRE(s.Kz.extent(0) == 60u);
  REQUIRE(s.qz.extent(0) == 60u);
  REQUIRE(s.h_pred.extent(0) == 60u);
}

TEST_CASE("deep copy isolates the copy from the original") {
  Grid g(4, 4, 1, 1.0, 1.0, 1.0);
  State s(g);
  Kokkos::deep_copy(s.eta, 3.0);

  RealArr1D copy("copy", g.nSurfaceCell());
  Kokkos::deep_copy(copy, s.eta);
  Kokkos::deep_copy(copy, 7.0);  // modify the copy only

  auto eta_h = Kokkos::create_mirror_view(s.eta);
  Kokkos::deep_copy(eta_h, s.eta);
  REQUIRE(eta_h(0) == Approx(3.0).margin(1e-12));  // original unchanged
}

TEST_CASE("time-level swap hn <- h, h <- h_new via deep_copy") {
  Grid g(2, 2, 3, 1.0, 1.0, 0.1, 1.0);
  GwState s(g);
  Kokkos::deep_copy(s.h, 2.0);
  Kokkos::deep_copy(s.hn, 0.0);

  RealArr1D h_new("h_new", g.nActiveCell());
  Kokkos::deep_copy(h_new, 5.0);

  Kokkos::deep_copy(s.hn, s.h);   // hn <- h
  Kokkos::deep_copy(s.h, h_new);  // h  <- h_new

  auto hn_h = Kokkos::create_mirror_view(s.hn);
  auto h_h = Kokkos::create_mirror_view(s.h);
  Kokkos::deep_copy(hn_h, s.hn);
  Kokkos::deep_copy(h_h, s.h);
  REQUIRE(hn_h(0) == Approx(2.0).margin(1e-12));
  REQUIRE(h_h(0) == Approx(5.0).margin(1e-12));
}
