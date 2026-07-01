#include "io/Hdf5Reader.hpp"

#include <hdf5.h>

#include <stdexcept>

#include "io/Hdf5Support.hpp"

namespace frehg2 {

namespace {

RealArr1DHost toView(const std::string& label, const std::vector<double>& v) {
  RealArr1DHost out(label, v.size());
  for (size_t i = 0; i < v.size(); ++i) out(i) = static_cast<real>(v[i]);
  return out;
}

h5::Guard openRead(const std::string& path) {
  if (H5Fis_hdf5(path.c_str()) <= 0) {
    throw std::runtime_error("Hdf5Reader: not a valid HDF5 file (or missing): '" + path +
                             "'");
  }
  hid_t f = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  h5::check(f, ("H5Fopen " + path).c_str());
  return h5::Guard(f, H5Fclose);
}

RealArr1DHost readField(const std::string& path, const std::string& group,
                        const std::string& field, long long time_step) {
  h5::Guard f = openRead(path);
  const std::string dset = group + "/" + field + "/" + std::to_string(time_step);
  if (!h5::linkExists(f.get(), dset)) {
    throw std::runtime_error("Hdf5Reader: dataset '" + dset + "' not found in '" + path +
                             "'");
  }
  return toView(field, h5::readDoubleDataset(f.get(), dset));
}

}  // namespace

RealArr1DHost Hdf5Reader::readSurface(const std::string& path, const std::string& field,
                                     long long time_step) {
  return readField(path, "surface", field, time_step);
}

RealArr1DHost Hdf5Reader::readSubsurface(const std::string& path, const std::string& field,
                                        long long time_step) {
  return readField(path, "subsurface", field, time_step);
}

RealArr1DHost Hdf5Reader::readSurfaceFilePerRank(const std::string& stem,
                                                 const std::string& field,
                                                 long long time_step, int nranks,
                                                 int global_size) {
  std::vector<double> global(static_cast<size_t>(global_size), 0.0);
  for (int r = 0; r < nranks; ++r) {
    const std::string path = stem + ".rank" + std::to_string(r) + ".h5";
    h5::Guard f = openRead(path);
    const std::string base = "surface/" + field + "/" + std::to_string(time_step);
    const std::string idxset = base + "_index";
    if (!h5::linkExists(f.get(), base)) continue;
    std::vector<double> vals = h5::readDoubleDataset(f.get(), base);
    std::vector<int> idx = h5::readIntDataset(f.get(), idxset);
    if (vals.size() != idx.size()) {
      throw std::runtime_error("Hdf5Reader: file_per_rank value/index size mismatch in " +
                               path);
    }
    for (size_t e = 0; e < vals.size(); ++e) {
      const int gi = idx[e];
      if (gi < 0 || gi >= global_size) {
        throw std::runtime_error("Hdf5Reader: file_per_rank index out of range in " + path);
      }
      global[static_cast<size_t>(gi)] = vals[e];
    }
  }
  return toView(field, global);
}

std::string Hdf5Reader::readMetadataAttr(const std::string& path, const std::string& attr) {
  h5::Guard f = openRead(path);
  h5::Guard g(H5Gopen2(f.get(), "simulation/metadata", H5P_DEFAULT), H5Gclose);
  h5::check(g.get(), "H5Gopen simulation/metadata");
  return h5::readStringAttr(g.get(), attr);
}

std::string Hdf5Reader::readSurfaceUnits(const std::string& path, const std::string& field,
                                        long long time_step) {
  h5::Guard f = openRead(path);
  const std::string dset = "surface/" + field + "/" + std::to_string(time_step);
  h5::Guard d(H5Dopen2(f.get(), dset.c_str(), H5P_DEFAULT), H5Dclose);
  h5::check(d.get(), ("H5Dopen " + dset).c_str());
  if (!h5::attrExists(d.get(), "units")) return "";
  return h5::readStringAttr(d.get(), "units");
}

RestartState Hdf5Reader::readCheckpoint(const std::string& path,
                                        const std::string& expected_config_sha256) {
  h5::Guard f = openRead(path);
  h5::Guard root(H5Gopen2(f.get(), "checkpoint", H5P_DEFAULT), H5Gclose);
  h5::check(root.get(), "H5Gopen checkpoint");

  RestartState rs;
  rs.time = h5::readDoubleAttr(root.get(), "time");
  rs.dt = h5::readDoubleAttr(root.get(), "dt");
  rs.storage_layout = h5::readStringAttr(root.get(), "storage_layout");
  rs.config_sha256 = h5::readStringAttr(root.get(), "config_sha256");
  rs.git_sha = h5::readStringAttr(root.get(), "git_sha");
  rs.has_sw = h5::readIntAttr(root.get(), "has_sw") != 0;
  rs.has_gw = h5::readIntAttr(root.get(), "has_gw") != 0;
  rs.has_solute = h5::readIntAttr(root.get(), "has_solute") != 0;

  auto rd = [&](const std::string& dset) {
    return toView("ckpt", h5::readDoubleDataset(f.get(), dset));
  };
  if (rs.has_sw) {
    rs.eta = rd("checkpoint/sw/eta");
    rs.u = rd("checkpoint/sw/u");
    rs.v = rd("checkpoint/sw/v");
  }
  if (rs.has_gw) {
    rs.h = rd("checkpoint/gw/h");
    rs.hn = rd("checkpoint/gw/hn");
    rs.wc = rd("checkpoint/gw/wc");
    rs.wcn = rd("checkpoint/gw/wcn");
  }
  if (rs.has_solute) {
    rs.C = rd("checkpoint/solute/C");
  }

  // Full solver state (P7.4): load every dataset under checkpoint/extra (if present).
  if (h5::linkExists(f.get(), "checkpoint/extra")) {
    h5::Guard grp(H5Gopen2(f.get(), "checkpoint/extra", H5P_DEFAULT), H5Gclose);
    h5::check(grp.get(), "H5Gopen checkpoint/extra");
    H5G_info_t info;
    h5::checkErr(H5Gget_info(grp.get(), &info), "H5Gget_info checkpoint/extra");
    for (hsize_t k = 0; k < info.nlinks; ++k) {
      ssize_t len = H5Lget_name_by_idx(grp.get(), ".", H5_INDEX_NAME, H5_ITER_INC, k, nullptr,
                                       0, H5P_DEFAULT);
      if (len < 0) throw std::runtime_error("Hdf5Reader: H5Lget_name_by_idx failed");
      std::string name(static_cast<size_t>(len), '\0');
      H5Lget_name_by_idx(grp.get(), ".", H5_INDEX_NAME, H5_ITER_INC, k, name.data(),
                         static_cast<size_t>(len) + 1, H5P_DEFAULT);
      rs.extra[name] = toView(name, h5::readDoubleDataset(f.get(), "checkpoint/extra/" + name));
    }
  }

  if (!expected_config_sha256.empty()) {
    rs.config_matches = (rs.config_sha256 == expected_config_sha256);
  }
  return rs;
}

}  // namespace frehg2
