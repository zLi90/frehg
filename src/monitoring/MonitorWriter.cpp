#include "monitoring/MonitorWriter.hpp"

#include <mpi.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>
#if !defined(_WIN32)
#include <unistd.h>
#endif

#include "monitoring/LineFlux.hpp"
#include "monitoring/ProbeLocator.hpp"
#include "re/GwFields.hpp"
#include "re/ReSolver.hpp"
#include "swe/SweFields.hpp"
#include "swe/SweSolver.hpp"

namespace frehg2 {

namespace {

void fsyncFile(FILE* fp) {
  if (fp == nullptr) return;
  std::fflush(fp);
#if !defined(_WIN32)
  std::fflush(fp);
  const int fd = fileno(fp);
  if (fd >= 0) fsync(fd);
#endif
}

}  // namespace

void MonitorWriter::configure(MonitorBundle bundle, std::string output_dir, std::string run_name,
                              int gnx, int gny, int gnz, double dx, double dy, double dz,
                              double x0, double y0, double botz, const MpiComm* mc) {
  bundle_ = std::move(bundle);
  output_dir_ = std::move(output_dir);
  run_name_ = std::move(run_name);
  gnx_ = gnx;
  gny_ = gny;
  gnz_ = gnz;
  dx_ = dx;
  dy_ = dy;
  dz_ = dz;
  x0_ = x0;
  y0_ = y0;
  botz_ = botz;
  mc_ = mc;

  buildProbeLocations(bundle_, gnx_, gny_, gnz_, dx_, dy_, dz_, x0_, y0_, botz_, mc_);

  columns_.clear();
  for (size_t pi = 0; pi < bundle_.probes.size(); ++pi) {
    for (const std::string& field : bundle_.probes[pi].fields) {
      Column col;
      col.kind = Column::Kind::Probe;
      col.probe_index = pi;
      col.header = bundle_.probes[pi].name + "." + field;
      columns_.push_back(std::move(col));
    }
  }
  for (size_t li = 0; li < bundle_.lines.size(); ++li) {
    Column col;
    col.kind = Column::Kind::Line;
    col.line_index = li;
    col.header = bundle_.lines[li].name + ".flux";
    columns_.push_back(std::move(col));
  }

  if (output_dir_.empty()) output_dir_ = ".";
  path_ = output_dir_ + "/monitors/" + run_name_ + ".csv";
}

void MonitorWriter::open(bool resume) {
  if (columns_.empty()) return;

  if (mc_ == nullptr || mc_->rank() == 0) {
    const std::filesystem::path dir = std::filesystem::path(path_).parent_path();
    if (!dir.empty()) std::filesystem::create_directories(dir);

    if (resume && std::filesystem::exists(path_)) {
      std::ifstream in(path_);
      std::string header;
      if (!std::getline(in, header))
        throw std::runtime_error("MonitorWriter: cannot read existing monitor CSV '" + path_ +
                                 "'");
      std::string expected = "time";
      for (const Column& c : columns_) expected += "," + c.header;
      if (header != expected)
        throw std::runtime_error("MonitorWriter: monitor CSV header mismatch on resume for '" +
                                 path_ + "'");
      header_written_ = true;
    } else {
      std::ofstream out(path_, std::ios::trunc);
      if (!out) throw std::runtime_error("MonitorWriter: cannot open '" + path_ + "' for write");
      out << "time";
      for (const Column& c : columns_) out << "," << c.header;
      out << "\n";
      out.flush();
      header_written_ = true;
    }
  }
  open_ = true;
}

void MonitorWriter::close() { open_ = false; }

double MonitorWriter::sampleProbeField(const ProbeSpec& probe, const std::string& field,
                                       const SweSolver* swe, const ReSolver* re,
                                       const Grid& swe_grid, const Grid& gw_grid,
                                       const RealArr1DHost* surf_conc,
                                       const RealArr1DHost* subs_conc) const {
  if (!probe.subsurface) {
    if (field == "C" || field == "conc" || field == "concentration") {
      if (surf_conc == nullptr || surf_conc->extent(0) == 0) return 0.0;
      const int c = swe_grid.getSurfaceIndex(probe.local_i, probe.local_j);
      return (*surf_conc)(static_cast<size_t>(c));
    }
    if (swe == nullptr) return 0.0;
    const int c = swe_grid.getSurfaceIndex(probe.local_i, probe.local_j);
    const SweFields& sf = swe->fields();
    if (field == "eta") return sf.eta(c);
    if (field == "u" || field == "uu") return sf.uu(c);
    if (field == "v" || field == "vv") return sf.vv(c);
    if (field == "dept" || field == "depth") return sf.dept(c);
    if (field == "bottom") return sf.bottom(c);
    throw std::runtime_error("MonitorWriter: unknown surface probe field '" + field + "'");
  }

  if (field == "C" || field == "conc" || field == "concentration") {
    if (subs_conc == nullptr || subs_conc->extent(0) == 0) return 0.0;
    const int c = gw_grid.getIndex(probe.local_i, probe.local_j, probe.local_k);
    return (*subs_conc)(static_cast<size_t>(c));
  }
  if (re == nullptr) return 0.0;
  const int c = gw_grid.getIndex(probe.local_i, probe.local_j, probe.local_k);
  const GwFields& gf = re->fields();
  if (field == "head" || field == "h") return gf.h(c);
  if (field == "moisture" || field == "wc") return gf.wc(c);
  if (field == "qx") return gf.qx(c);
  if (field == "qy") return gf.qy(c);
  if (field == "qz") return gf.qz(c);
  throw std::runtime_error("MonitorWriter: unknown subsurface probe field '" + field + "'");
}

void MonitorWriter::writeRow(real time, const SweSolver* swe, const ReSolver* re,
                             const Grid& swe_grid, const Grid& gw_grid,
                             const RealArr1DHost* surf_conc, const RealArr1DHost* subs_conc) {
  if (!open_ || columns_.empty()) return;

  const int rank = mc_ ? mc_->rank() : 0;
  const int size = mc_ ? mc_->size() : 1;
  const int ncols = static_cast<int>(columns_.size());
  std::vector<double> local_vals(static_cast<size_t>(ncols), 0.0);

  for (int ci = 0; ci < ncols; ++ci) {
    const Column& col = columns_[static_cast<size_t>(ci)];
    if (col.kind == Column::Kind::Line) {
      const LineFluxSpec& line = bundle_.lines[col.line_index];
      local_vals[static_cast<size_t>(ci)] = integrateLineFlux(
          line, swe, re, swe_grid, gw_grid, gnx_, gny_, dx_, dy_, x0_, y0_, mc_);
      if (size > 1) {
        double global = 0.0;
        MPI_Allreduce(&local_vals[static_cast<size_t>(ci)], &global, 1, MPI_DOUBLE, MPI_SUM,
                      mc_->comm());
        local_vals[static_cast<size_t>(ci)] = global;
      }
      continue;
    }

    const ProbeSpec& probe = bundle_.probes[col.probe_index];
    const std::string want = col.header.substr(probe.name.size() + 1);
    double val = 0.0;
    if (rank == probe.owner_rank) {
      val = sampleProbeField(probe, want, swe, re, swe_grid, gw_grid, surf_conc, subs_conc);
    }
    if (size > 1) {
      double global = 0.0;
      MPI_Allreduce(&val, &global, 1, MPI_DOUBLE, MPI_SUM, mc_->comm());
      val = global;
    }
    local_vals[static_cast<size_t>(ci)] = val;
  }

  if (rank != 0) return;

  FILE* fp = std::fopen(path_.c_str(), "a");
  if (fp == nullptr)
    throw std::runtime_error("MonitorWriter: cannot append to '" + path_ + "'");
  std::fprintf(fp, "%.9g", static_cast<double>(time));
  for (double v : local_vals) std::fprintf(fp, ",%.9g", v);
  std::fprintf(fp, "\n");
  fsyncFile(fp);
  std::fclose(fp);
}

}  // namespace frehg2
