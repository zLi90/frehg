// Shared b1-sw driver for serial (4c) and MPI rank-equivalence (4e) gates.
#pragma once

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "linear/DomainDecomposition.hpp"
#include "linear/backends/PetscLinearSolver.hpp"
#include "swe/SweSolver.hpp"

namespace frehg2 {
namespace b1_sw {

inline std::vector<double> readDoubles(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::vector<double> v;
  double x;
  while (in >> x) v.push_back(x);
  return v;
}

inline double interpBc(const std::vector<double>& t, const std::vector<double>& val,
                       double tc) {
  const int n = static_cast<int>(t.size());
  if (tc == 0.0) return val[0];
  int ind = 1;
  if (tc > t[ind]) {
    ind += 1;
    while (ind < n && tc > t[ind]) ind += 1;
  }
  if (ind >= n) return val[n - 1];
  return val[ind - 1] + (val[ind] - val[ind - 1]) * (tc - t[ind - 1]) / (t[ind] - t[ind - 1]);
}

inline double relL2(const std::vector<double>& a, const std::vector<double>& b) {
  double num = 0.0, den = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    const double d = a[i] - b[i];
    num += d * d;
    den += b[i] * b[i];
  }
  if (den == 0.0) return std::sqrt(num);
  return std::sqrt(num) / std::sqrt(den);
}

struct B1Params {
  int gnx = 1;
  int gny = 10;
  double dx = 80.0;
  double dy = 80.0;
  double dt = 5.0;
  int nsteps = 3600;
  int out_every = 360;
  std::vector<double> bath;
  std::vector<double> rain_t;
  std::vector<double> rain_v;
  double offset = 0.0;
};

inline B1Params loadB1Params(const std::string& legacy_dir) {
  B1Params p;
  p.bath = readDoubles(legacy_dir + "/b1-input/bath");
  p.gny = static_cast<int>(p.bath.size());
  const auto rainflat = readDoubles(legacy_dir + "/b1-input/rain");
  if (rainflat.size() != 8u) throw std::runtime_error("unexpected rain file size");
  p.rain_t.resize(4);
  p.rain_v.resize(4);
  for (int k = 0; k < 4; ++k) {
    p.rain_t[static_cast<size_t>(k)] = rainflat[2 * k];
    p.rain_v[static_cast<size_t>(k)] = rainflat[2 * k + 1];
  }
  double bmin = p.bath[0];
  for (double b : p.bath) bmin = std::min(bmin, b);
  p.offset = -bmin;
  return p;
}

struct B1Snapshot {
  int t = 0;
  std::vector<double> depth;
};

struct B1RunResult {
  std::vector<B1Snapshot> snapshots;
  double aggregate_rel_l2_vs_legacy = 0.0;
  double worst_snapshot_rel_l2 = 0.0;
  int worst_t = 0;
};

// Run the full b1-sw case with the given MPI decomposition. Returns gathered global
// depth snapshots on rank 0 (empty on other ranks).
inline B1RunResult runB1(const B1Params& p, MpiComm& mc, const std::string& legacy_dir) {
  B1RunResult out;
  if (mc.rank() != 0) out.snapshots.clear();

  Grid grid(mc.localNx(), mc.localNy(), 1, p.dx, p.dy, 0.1);
  SweSolver swe(grid, &mc);
  SweParams sp;
  sp.gravity = 9.81;
  sp.manning = 0.019;
  sp.min_depth = 1.0e-8;
  sp.visc_x = 1.0e-6;
  sp.visc_y = 1.0e-6;
  sp.hD = 0.1;
  sp.wtfh = 1.0e-8;
  sp.dt = p.dt;
  sp.offset = p.offset;
  swe.setParams(sp);

  RealArr1DHost bed_global("bed_global",
                           static_cast<size_t>(p.gnx * p.gny));
  for (int j = 0; j < p.gny; ++j)
    for (int i = 0; i < p.gnx; ++i)
      bed_global(static_cast<size_t>(i + j * p.gnx)) =
          p.bath[static_cast<size_t>(j)] + p.offset;
  swe.setBathymetry(bed_global);
  swe.initializeState(-2.0);

  Decomp2D dd(mc);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1e-12;
  PetscLinearSolver solver(cfg);
  swe.attachSolver(solver, dd);

  RealArr1DHost global_dept("global_dept",
                            static_cast<size_t>(p.gnx * p.gny));
  double agg_num = 0.0, agg_den = 0.0;
  double worst = 0.0;
  int worst_t = 0;
  double t_current = 0.0;

  for (int step = 1; step <= p.nsteps; ++step) {
    t_current += p.dt;
    const double rain = interpBc(p.rain_t, p.rain_v, t_current);
    swe.advanceStep(rain, 0.0);

    if (step % p.out_every == 0) {
      const int tsave = static_cast<int>(std::llround(t_current));
      swe.gatherDeptGlobal(global_dept);
      if (mc.rank() == 0) {
        std::vector<double> ref =
            readDoubles(legacy_dir + "/reference/depth_" + std::to_string(tsave));
        std::vector<double> got(static_cast<size_t>(p.gnx * p.gny));
        for (int j = 0; j < p.gny; ++j)
          for (int i = 0; i < p.gnx; ++i)
            got[static_cast<size_t>(i + j * p.gnx)] =
                global_dept(static_cast<size_t>(i + j * p.gnx));
        const double e = relL2(got, ref);
        for (size_t k = 0; k < got.size(); ++k) {
          const double d = got[k] - ref[k];
          agg_num += d * d;
          agg_den += ref[k] * ref[k];
        }
        if (e > worst) {
          worst = e;
          worst_t = tsave;
        }
        out.snapshots.push_back({tsave, got});
      }
    }
  }

  if (mc.rank() == 0) {
    out.aggregate_rel_l2_vs_legacy = std::sqrt(agg_num) / std::sqrt(agg_den);
    out.worst_snapshot_rel_l2 = worst;
    out.worst_t = worst_t;
  }
  return out;
}

inline void procGridY(int size, int& mpi_nx, int& mpi_ny) {
  // b1-sw is NX=1; decompose only in y.
  mpi_nx = 1;
  mpi_ny = size;
}

}  // namespace b1_sw
}  // namespace frehg2
