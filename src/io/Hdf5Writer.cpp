#include "io/Hdf5Writer.hpp"

#include <sys/stat.h>  // mkdir
#include <unistd.h>    // fsync

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "core/Grid.hpp"
#include "io/Hdf5Support.hpp"

namespace frehg2 {

namespace {

// Dataset/checkpoint key for a snapshot time. Whole-second times keep their plain integer
// name ("0", "1800") for backward compatibility; fractional (sub-second) times get a lossless,
// trimmed decimal name ("0.3", "0.0001") so two distinct snapshots never collide.
std::string timeKey(double t) {
  const long long r = std::llround(t);
  if (std::fabs(t - static_cast<double>(r)) < 1.0e-9) {
    return std::to_string(r);
  }
  std::ostringstream os;
  os << std::setprecision(12) << t;
  return os.str();
}

std::string stemOf(const std::string& path) {
  const auto dot = path.find_last_of('.');
  const auto slash = path.find_last_of('/');
  if (dot != std::string::npos && (slash == std::string::npos || dot > slash)) {
    return path.substr(0, dot);
  }
  return path;
}

std::string dirOf(const std::string& path) {
  const auto slash = path.find_last_of('/');
  return slash == std::string::npos ? std::string(".") : path.substr(0, slash);
}

std::string isoNow() {
  std::time_t t = std::time(nullptr);
  char buf[32];
  std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", std::localtime(&t));
  return std::string(buf);
}

std::string defaultUnits(const std::string& name) {
  if (name == "water_depth" || name == "water_surface_elevation" || name == "eta" ||
      name == "hydraulic_head") {
    return "m";
  }
  if (name == "velocity_x" || name == "velocity_y" || name == "uu" || name == "vv") {
    return "m/s";
  }
  if (name == "water_content") return "m^3/m^3";
  if (name.rfind("darcy_flux", 0) == 0) return "m/s";
  return "";
}

}  // namespace

Hdf5Writer::Hdf5Writer(MPI_Comm comm, IoLayout layout, IoMode mode, int max_checkpoints)
    : comm_(comm), layout_(std::move(layout)), mode_(mode), max_checkpoints_(max_checkpoints) {
  MPI_Comm_rank(comm_, &rank_);
  MPI_Comm_size(comm_, &size_);
}

Hdf5Writer::~Hdf5Writer() {
  if (file_ >= 0) {
    H5Fclose(file_);
    file_ = -1;
  }
}

std::string Hdf5Writer::rankFilePath() const {
  return file_stem_ + ".rank" + std::to_string(rank_) + ".h5";
}

void Hdf5Writer::openFile(const std::string& path, const RunMetadata& meta) {
  path_ = path;
  file_stem_ = stemOf(path);
  meta_ = meta;
  meta_.mpi_ranks = size_;
  if (meta_.date.empty()) meta_.date = isoNow();

  // Make sure the output directory exists (best-effort; rank 0 only).
  if (rank_ == 0) {
    const std::string dir = dirOf(path);
    if (dir != ".") {
      // POSIX mkdir -p semantics via a tiny loop.
      std::string acc;
      for (size_t i = 0; i < dir.size(); ++i) {
        acc.push_back(dir[i]);
        if (dir[i] == '/' || i + 1 == dir.size()) {
          if (!acc.empty() && acc != "/") {
            ::mkdir(acc.c_str(), 0755);
          }
        }
      }
    }
  }
  MPI_Barrier(comm_);

  switch (mode_) {
    case IoMode::SerialGather: {
      if (rank_ == 0) {
        file_ = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        h5::check(file_, "H5Fcreate (serial_gather)");
        owns_file_ = true;
        writeMetadata(file_);
      }
      break;
    }
    case IoMode::ParallelCollective: {
      hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
      h5::check(fapl, "H5Pcreate fapl");
      h5::checkErr(H5Pset_fapl_mpio(fapl, comm_, MPI_INFO_NULL), "H5Pset_fapl_mpio");
      file_ = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
      H5Pclose(fapl);
      h5::check(file_, "H5Fcreate (parallel_collective)");
      owns_file_ = true;
      writeMetadata(file_);  // collective; all ranks write identical attributes
      break;
    }
    case IoMode::FilePerRank: {
      file_ = H5Fcreate(rankFilePath().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      h5::check(file_, "H5Fcreate (file_per_rank)");
      owns_file_ = true;
      writeMetadata(file_);
      break;
    }
  }
}

void Hdf5Writer::writeMetadata(hid_t file) {
  h5::Guard meta = h5::ensureGroup(file, "simulation/metadata");
  hid_t g = meta.get();
  h5::writeStringAttr(g, "title", meta_.title);
  h5::writeStringAttr(g, "version", meta_.version);
  h5::writeStringAttr(g, "date", meta_.date);
  h5::writeStringAttr(g, "frehg2_version", meta_.frehg2_version);
  h5::writeStringAttr(g, "git_sha", meta_.git_sha);
  h5::writeIntAttr(g, "git_dirty", meta_.git_dirty ? 1 : 0);
  h5::writeStringAttr(g, "config_sha256", meta_.config_sha256);
  h5::writeStringAttr(g, "build_type", meta_.build_type);
  h5::writeStringAttr(g, "compiler", meta_.compiler);
  h5::writeStringAttr(g, "kokkos_backend", meta_.kokkos_backend);
  h5::writeStringAttr(g, "solver_backend", meta_.solver_backend);
  h5::writeIntAttr(g, "mpi_ranks", meta_.mpi_ranks);
  h5::writeStringAttr(g, "io_mode", ioModeName(mode_));
}

void Hdf5Writer::writeDomain(const Grid& grid) {
  grid_nx_ = grid.nx();
  grid_ny_ = grid.ny();
  grid_nz_ = grid.nz();
  grid_written_ = true;

  const bool do_write = (mode_ == IoMode::ParallelCollective) || owns_file_;
  if (!do_write || file_ < 0) return;

  h5::Guard dom = h5::ensureGroup(file_, "domain");
  hid_t g = dom.get();
  h5::writeIntAttr(g, "nx", grid.nx());
  h5::writeIntAttr(g, "ny", grid.ny());
  h5::writeIntAttr(g, "nz", grid.nz());
  h5::writeDoubleAttr(g, "dx", static_cast<double>(grid.dx()));
  h5::writeDoubleAttr(g, "dy", static_cast<double>(grid.dy()));
  h5::writeDoubleAttr(g, "dz", static_cast<double>(grid.dz()));
  h5::writeDoubleAttr(g, "origin_x", 0.0);
  h5::writeDoubleAttr(g, "origin_y", 0.0);
  h5::writeDoubleAttr(g, "origin_z", 0.0);

  const double dims[3] = {static_cast<double>(grid.nx()), static_cast<double>(grid.ny()),
                          static_cast<double>(grid.nz())};
  h5::writeDoubleDataset(g, "dimensions", dims, 3, false);
}

void Hdf5Writer::writeArrayWithUnits(hid_t loc, const std::string& dsname, const double* data,
                                     hsize_t n, const std::string& units, bool gzip) {
  hsize_t dims[1] = {n};
  h5::Guard space(H5Screate_simple(1, dims, nullptr), H5Sclose);
  h5::Guard dcpl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose);
  if (gzip && n > 0) {
    hsize_t chunk[1] = {n};
    h5::checkErr(H5Pset_chunk(dcpl.get(), 1, chunk), "H5Pset_chunk");
    h5::checkErr(H5Pset_deflate(dcpl.get(), 6), "H5Pset_deflate");
  }
  h5::Guard dset(H5Dcreate2(loc, dsname.c_str(), H5T_IEEE_F64LE, space.get(), H5P_DEFAULT,
                            dcpl.get(), H5P_DEFAULT),
                 H5Dclose);
  h5::check(dset.get(), ("H5Dcreate " + dsname).c_str());
  if (n > 0) {
    h5::checkErr(H5Dwrite(dset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data),
                 ("H5Dwrite " + dsname).c_str());
  }
  if (!units.empty()) h5::writeStringAttr(dset.get(), "units", units);
}

std::vector<double> Hdf5Writer::gatherGlobal(const RealArr1DHost& owned,
                                             const IntArr1DHost& global_idx,
                                             int global_size) const {
  const int n_local = static_cast<int>(owned.extent(0));
  std::vector<int> counts(static_cast<size_t>(size_), 0);
  MPI_Gather(&n_local, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm_);

  std::vector<int> displs(static_cast<size_t>(size_), 0);
  int total = 0;
  if (rank_ == 0) {
    for (int r = 0; r < size_; ++r) {
      displs[static_cast<size_t>(r)] = total;
      total += counts[static_cast<size_t>(r)];
    }
  }

  std::vector<double> local_vals(static_cast<size_t>(n_local));
  std::vector<int> local_idx(static_cast<size_t>(n_local));
  for (int i = 0; i < n_local; ++i) {
    local_vals[static_cast<size_t>(i)] = static_cast<double>(owned(static_cast<size_t>(i)));
    local_idx[static_cast<size_t>(i)] = global_idx(static_cast<size_t>(i));
  }

  std::vector<double> all_vals;
  std::vector<int> all_idx;
  if (rank_ == 0) {
    all_vals.resize(static_cast<size_t>(total));
    all_idx.resize(static_cast<size_t>(total));
  }
  MPI_Gatherv(local_vals.data(), n_local, MPI_DOUBLE,
              rank_ == 0 ? all_vals.data() : nullptr, counts.data(), displs.data(),
              MPI_DOUBLE, 0, comm_);
  MPI_Gatherv(local_idx.data(), n_local, MPI_INT, rank_ == 0 ? all_idx.data() : nullptr,
              counts.data(), displs.data(), MPI_INT, 0, comm_);

  std::vector<double> global;
  if (rank_ == 0) {
    global.assign(static_cast<size_t>(global_size), 0.0);
    for (int e = 0; e < total; ++e) {
      const int gi = all_idx[static_cast<size_t>(e)];
      if (gi < 0 || gi >= global_size) {
        throw std::runtime_error("Hdf5Writer: global index out of range during gather");
      }
      global[static_cast<size_t>(gi)] = all_vals[static_cast<size_t>(e)];
    }
  }
  return global;
}

void Hdf5Writer::writeFieldSerialGather(const std::string& group, const std::string& name,
                                        double time, const RealArr1DHost& owned,
                                        const IntArr1DHost& gidx, int global_size,
                                        const std::string& units) {
  std::vector<double> global = gatherGlobal(owned, gidx, global_size);
  if (rank_ != 0) return;
  h5::Guard grp = h5::ensureGroup(file_, group + "/" + name);
  writeArrayWithUnits(grp.get(), timeKey(time), global.data(),
                      static_cast<hsize_t>(global_size), units, /*gzip=*/true);
}

void Hdf5Writer::writeFieldParallel(const std::string& group, const std::string& name,
                                    double time, const RealArr1DHost& owned,
                                    const IntArr1DHost& gidx, int global_size,
                                    const std::string& units) {
  // Collective: every rank creates the dataset and writes its owned elements by global
  // index. Compression disabled (collective filtered writes are avoided per P3.2.1a).
  h5::Guard grp = h5::ensureGroup(file_, group + "/" + name);
  const hsize_t gdim = static_cast<hsize_t>(global_size);
  h5::Guard fspace(H5Screate_simple(1, &gdim, nullptr), H5Sclose);
  h5::Guard dset(H5Dcreate2(grp.get(), timeKey(time).c_str(), H5T_IEEE_F64LE, fspace.get(),
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                 H5Dclose);
  h5::check(dset.get(), "H5Dcreate (parallel)");

  const int n_local = static_cast<int>(owned.extent(0));
  std::vector<double> vals(static_cast<size_t>(n_local));
  std::vector<hsize_t> coords(static_cast<size_t>(n_local));
  for (int i = 0; i < n_local; ++i) {
    vals[static_cast<size_t>(i)] = static_cast<double>(owned(static_cast<size_t>(i)));
    coords[static_cast<size_t>(i)] = static_cast<hsize_t>(gidx(static_cast<size_t>(i)));
  }

  h5::Guard filespace(H5Dget_space(dset.get()), H5Sclose);
  if (n_local > 0) {
    h5::checkErr(H5Sselect_elements(filespace.get(), H5S_SELECT_SET,
                                    static_cast<size_t>(n_local), coords.data()),
                 "H5Sselect_elements");
  } else {
    h5::checkErr(H5Sselect_none(filespace.get()), "H5Sselect_none");
  }
  const hsize_t mdim = static_cast<hsize_t>(n_local);
  h5::Guard mspace(H5Screate_simple(1, &mdim, nullptr), H5Sclose);
  if (n_local == 0) h5::checkErr(H5Sselect_none(mspace.get()), "H5Sselect_none mem");

  h5::Guard dxpl(H5Pcreate(H5P_DATASET_XFER), H5Pclose);
  h5::checkErr(H5Pset_dxpl_mpio(dxpl.get(), H5FD_MPIO_COLLECTIVE), "H5Pset_dxpl_mpio");
  h5::checkErr(H5Dwrite(dset.get(), H5T_NATIVE_DOUBLE, mspace.get(), filespace.get(),
                        dxpl.get(), vals.data()),
               "H5Dwrite (parallel)");
  if (!units.empty()) h5::writeStringAttr(dset.get(), "units", units);
}

void Hdf5Writer::writeFieldPerRank(const std::string& group, const std::string& name,
                                   double time, const RealArr1DHost& owned,
                                   const IntArr1DHost& gidx, const std::string& units) {
  // Each rank writes its owned values + the global-index map so a reader can reassemble.
  h5::Guard grp = h5::ensureGroup(file_, group + "/" + name);
  const int n_local = static_cast<int>(owned.extent(0));
  std::vector<double> vals(static_cast<size_t>(n_local));
  std::vector<int> idx(static_cast<size_t>(n_local));
  for (int i = 0; i < n_local; ++i) {
    vals[static_cast<size_t>(i)] = static_cast<double>(owned(static_cast<size_t>(i)));
    idx[static_cast<size_t>(i)] = gidx(static_cast<size_t>(i));
  }
  writeArrayWithUnits(grp.get(), timeKey(time), vals.data(),
                      static_cast<hsize_t>(n_local), units, /*gzip=*/true);
  h5::writeIntDataset(grp.get(), timeKey(time) + "_index", idx.data(),
                      static_cast<hsize_t>(n_local));
}

void Hdf5Writer::writeSurfaceField(const std::string& name, double time,
                                   const RealArr1DHost& owned_values,
                                   const std::string& units) {
  const std::string u = units.empty() ? defaultUnits(name) : units;
  const int gsize = layout_.nx * layout_.ny;
  switch (mode_) {
    case IoMode::SerialGather:
      writeFieldSerialGather("surface", name, time, owned_values, layout_.surf_global_idx,
                             gsize, u);
      break;
    case IoMode::ParallelCollective:
      writeFieldParallel("surface", name, time, owned_values, layout_.surf_global_idx, gsize,
                         u);
      break;
    case IoMode::FilePerRank:
      writeFieldPerRank("surface", name, time, owned_values, layout_.surf_global_idx, u);
      break;
  }
  surf_times_[name][time] = timeKey(time);
}

void Hdf5Writer::writeSubsurfaceField(const std::string& name, double time,
                                      const RealArr1DHost& owned_values,
                                      const std::string& units) {
  const std::string u = units.empty() ? defaultUnits(name) : units;
  const int gsize = layout_.nx * layout_.ny * layout_.nz;
  switch (mode_) {
    case IoMode::SerialGather:
      writeFieldSerialGather("subsurface", name, time, owned_values,
                             layout_.subs_global_idx, gsize, u);
      break;
    case IoMode::ParallelCollective:
      writeFieldParallel("subsurface", name, time, owned_values, layout_.subs_global_idx,
                         gsize, u);
      break;
    case IoMode::FilePerRank:
      writeFieldPerRank("subsurface", name, time, owned_values, layout_.subs_global_idx, u);
      break;
  }
  subs_times_[name][time] = timeKey(time);
}

void Hdf5Writer::writeMonitorSeries(const std::string& name, const std::vector<double>& times,
                                    const std::vector<double>& values) {
  if (rank_ != 0) return;  // monitors are global series held by rank 0
  if (times.size() != values.size()) {
    throw std::runtime_error("Hdf5Writer::writeMonitorSeries: times/values size mismatch");
  }
  // Serial-mode files are owned by rank 0; for the shared parallel file, write monitors to
  // a separate rank-0 file to avoid collective-metadata constraints.
  const bool to_main = owns_file_ && mode_ != IoMode::ParallelCollective;
  hid_t target = -1;
  hid_t mon_file = -1;
  if (to_main) {
    target = file_;
  } else {
    mon_file = H5Fcreate((file_stem_ + ".monitor.h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                         H5P_DEFAULT);
    h5::check(mon_file, "H5Fcreate monitor");
    target = mon_file;
  }
  {
    h5::Guard grp = h5::ensureGroup(target, "monitoring/" + name);
    h5::writeDoubleDataset(grp.get(), "time", times.data(),
                           static_cast<hsize_t>(times.size()), false);
    h5::writeDoubleDataset(grp.get(), "value", values.data(),
                           static_cast<hsize_t>(values.size()), false);
  }
  if (mon_file >= 0) H5Fclose(mon_file);
}

void Hdf5Writer::writeCheckpoint(const CheckpointState& cp) {
  // Atomic, separate checkpoint file written by rank 0: temp file -> flush+fsync -> rename.
  if (rank_ != 0) {
    MPI_Barrier(comm_);
    return;
  }
  const std::string final_path =
      file_stem_ + ".ckpt." + timeKey(cp.time) + ".h5";
  const std::string tmp_path = final_path + ".tmp";

  {
    hid_t f = H5Fcreate(tmp_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    h5::check(f, "H5Fcreate checkpoint tmp");
    h5::Guard file(f, H5Fclose);

    h5::Guard root = h5::ensureGroup(f, "checkpoint");
    h5::writeDoubleAttr(root.get(), "time", cp.time);
    h5::writeDoubleAttr(root.get(), "dt", cp.dt);
    h5::writeStringAttr(root.get(), "storage_layout", cp.storage_layout);
    h5::writeStringAttr(root.get(), "config_sha256", cp.config_sha256);
    h5::writeStringAttr(root.get(), "git_sha", cp.git_sha);
    h5::writeIntAttr(root.get(), "has_sw", cp.has_sw ? 1 : 0);
    h5::writeIntAttr(root.get(), "has_gw", cp.has_gw ? 1 : 0);
    h5::writeIntAttr(root.get(), "has_solute", cp.has_solute ? 1 : 0);

    auto writeArr = [&](const std::string& grp, const std::string& nm,
                        const RealArr1DHost& a) {
      h5::Guard g = h5::ensureGroup(f, grp);
      std::vector<double> buf(a.extent(0));
      for (size_t i = 0; i < a.extent(0); ++i) buf[i] = static_cast<double>(a(i));
      h5::writeDoubleDataset(g.get(), nm, buf.data(), static_cast<hsize_t>(buf.size()),
                             false);
    };
    if (cp.has_sw) {
      writeArr("checkpoint/sw", "eta", cp.eta);
      writeArr("checkpoint/sw", "u", cp.u);
      writeArr("checkpoint/sw", "v", cp.v);
    }
    if (cp.has_gw) {
      writeArr("checkpoint/gw", "h", cp.h);
      writeArr("checkpoint/gw", "hn", cp.hn);
      writeArr("checkpoint/gw", "wc", cp.wc);
      writeArr("checkpoint/gw", "wcn", cp.wcn);
    }
    if (cp.has_solute) {
      writeArr("checkpoint/solute", "C", cp.C);
    }
    // Full solver state (P7.4): every halo-padded field, for bit-exact restart.
    for (const auto& kv : cp.extra) {
      writeArr("checkpoint/extra", kv.first, kv.second);
    }
    h5::checkErr(H5Fflush(f, H5F_SCOPE_GLOBAL), "H5Fflush checkpoint");
  }  // file closed here

  // fsync the temp file, then atomically rename over the final path.
  {
    FILE* fp = std::fopen(tmp_path.c_str(), "rb");
    if (fp) {
      ::fsync(fileno(fp));
      std::fclose(fp);
    }
  }
  if (std::rename(tmp_path.c_str(), final_path.c_str()) != 0) {
    throw std::runtime_error("Hdf5Writer: checkpoint rename failed for " + final_path);
  }

  // Prune to the most recent max_checkpoints_ files.
  checkpoint_times_.push_back(cp.time);
  if (static_cast<int>(checkpoint_times_.size()) > max_checkpoints_) {
    const double old = checkpoint_times_.front();
    checkpoint_times_.erase(checkpoint_times_.begin());
    const std::string old_path = file_stem_ + ".ckpt." + timeKey(old) + ".h5";
    std::remove(old_path.c_str());
  }
  MPI_Barrier(comm_);
}

void Hdf5Writer::writeXdmf() {
  if (rank_ != 0 || !grid_written_) return;
  // Reference the canonical file: serial_gather/parallel_collective share `path_`.
  std::string h5name = path_;
  if (mode_ == IoMode::FilePerRank) {
    // file_per_rank stores per-rank slabs; a VDS index would be needed for a single
    // sidecar. Reference rank 0's file as a smoke-level sidecar.
    h5name = rankFilePath();
  }
  const auto slash = h5name.find_last_of('/');
  const std::string h5rel = slash == std::string::npos ? h5name : h5name.substr(slash + 1);

  // Build a temporal collection from the surface fields we actually wrote (time -> key).
  std::map<double, std::string> all_times;
  for (const auto& kv : surf_times_) {
    for (const auto& tk : kv.second) all_times[tk.first] = tk.second;
  }

  std::ofstream xmf(file_stem_ + ".xmf");
  if (!xmf) return;
  const int nx = grid_nx_;
  const int ny = grid_ny_;
  xmf << "<?xml version=\"1.0\" ?>\n";
  xmf << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
  xmf << "<Xdmf Version=\"3.0\">\n  <Domain>\n";
  xmf << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" "
         "CollectionType=\"Temporal\">\n";
  for (const auto& tkey : all_times) {
    const double t = tkey.first;
    const std::string& key = tkey.second;
    xmf << "      <Grid Name=\"t" << key << "\" GridType=\"Uniform\">\n";
    xmf << "        <Time Value=\"" << t << "\"/>\n";
    xmf << "        <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << ny << " " << nx
        << "\"/>\n";
    xmf << "        <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
    xmf << "          <DataItem Dimensions=\"2\" Format=\"XML\">0 0</DataItem>\n";
    xmf << "          <DataItem Dimensions=\"2\" Format=\"XML\">1 1</DataItem>\n";
    xmf << "        </Geometry>\n";
    for (const auto& kv : surf_times_) {
      if (kv.second.count(t) == 0) continue;
      xmf << "        <Attribute Name=\"" << kv.first
          << "\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
      xmf << "          <DataItem Dimensions=\"" << ny << " " << nx
          << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << h5rel
          << ":/surface/" << kv.first << "/" << key << "</DataItem>\n";
      xmf << "        </Attribute>\n";
    }
    xmf << "      </Grid>\n";
  }
  xmf << "    </Grid>\n  </Domain>\n</Xdmf>\n";
}

void Hdf5Writer::close() {
  if (file_ >= 0) {
    H5Fflush(file_, H5F_SCOPE_GLOBAL);
  }
  writeXdmf();
  if (file_ >= 0) {
    H5Fclose(file_);
    file_ = -1;
  }
  owns_file_ = false;
}

}  // namespace frehg2
