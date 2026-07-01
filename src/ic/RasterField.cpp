#include "ic/RasterField.hpp"

#include <fstream>
#include <stdexcept>

#include "io/AsciiRaster.hpp"
#include "io/Config.hpp"
#include "io/Hdf5Support.hpp"

namespace frehg2 {

namespace {

std::vector<double> readDoublesList(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("RasterField: cannot open '" + path + "'");
  std::vector<double> raw;
  double v;
  while (in >> v) raw.push_back(v);
  return raw;
}

RasterField fromDoubles(const std::vector<double>& raw, int gnx, int gny, int gnz,
                        const std::string& path) {
  RasterField rf;
  rf.gnx = gnx;
  rf.gny = gny;
  rf.gnz = gnz;
  const size_t expect2 = static_cast<size_t>(gnx) * static_cast<size_t>(gny);
  const size_t expect3 = expect2 * static_cast<size_t>(gnz);
  if (raw.size() == expect3) {
    rf.values = raw;
    return rf;
  }
  if (raw.size() == expect2) {
    rf.gnz = 1;
    rf.values = raw;
    return rf;
  }
  if (gnx == 1 && raw.size() == static_cast<size_t>(gny)) {
    rf.values.resize(expect2);
    for (int j = 0; j < gny; ++j) rf.values[static_cast<size_t>(j)] = raw[static_cast<size_t>(j)];
    rf.gnz = 1;
    return rf;
  }
  throw std::runtime_error("RasterField: '" + path + "' has " + std::to_string(raw.size()) +
                         " values; expected " + std::to_string(gnx * gny) + " or " +
                         std::to_string(gnx * gny * gnz));
}

RasterField loadHdf5(const std::string& path, const std::string& dataset, int gnx, int gny,
                     int gnz) {
  h5::Guard file(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
  if (!file.valid()) throw std::runtime_error("RasterField: cannot open HDF5 '" + path + "'");
  h5::Guard dset(H5Dopen2(file.get(), dataset.c_str(), H5P_DEFAULT), H5Dclose);
  if (!dset.valid())
    throw std::runtime_error("RasterField: dataset '" + dataset + "' not found in '" + path + "'");
  h5::Guard space(H5Dget_space(dset.get()), H5Sclose);
  const int rank = H5Sget_simple_extent_ndims(space.get());
  if (rank < 1 || rank > 3)
    throw std::runtime_error("RasterField: dataset '" + dataset + "' must be rank 1–3");
  std::vector<hsize_t> dims(static_cast<size_t>(rank));
  H5Sget_simple_extent_dims(space.get(), dims.data(), nullptr);
  size_t n = 1;
  for (hsize_t d : dims) n *= static_cast<size_t>(d);
  std::vector<double> raw(n);
  h5::checkErr(H5Dread(dset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw.data()),
               "H5Dread");
  // Accept (nx,ny), (ny,nx) with ny==1, or (nx,ny,nz).
  if (rank == 2 && static_cast<int>(dims[0]) == gnx && static_cast<int>(dims[1]) == gny)
    return fromDoubles(raw, gnx, gny, gnz, path);
  if (rank == 3 && static_cast<int>(dims[0]) == gnx && static_cast<int>(dims[1]) == gny &&
      static_cast<int>(dims[2]) == gnz) {
    RasterField rf;
    rf.gnx = gnx;
    rf.gny = gny;
    rf.gnz = gnz;
    rf.values = std::move(raw);
    return rf;
  }
  if (rank == 1 && static_cast<int>(n) == gnx * gny * gnz) return fromDoubles(raw, gnx, gny, gnz, path);
  throw std::runtime_error("RasterField: HDF5 '" + dataset + "' shape does not match domain");
}

}  // namespace

RasterField loadRasterField(const Config& config, const std::string& file,
                            const std::string& format, const std::string& dataset, int gnx,
                            int gny, int gnz) {
  const std::string path = config.resolvePath(file);
  std::string fmt = format;
  bool as_raster = (fmt == "raster" || fmt == "esri" || fmt == "ascii_raster");
  bool as_hdf5 = (fmt == "hdf5");
  if (fmt == "auto") {
    if (path.size() >= 3 && path.compare(path.size() - 3, 3, ".h5") == 0) as_hdf5 = true;
    else if (path.size() >= 4 && path.compare(path.size() - 4, 4, ".asc") == 0) as_raster = true;
  }
  if (as_hdf5) return loadHdf5(path, dataset, gnx, gny, gnz);
  if (as_raster) {
    const AsciiRaster ras = AsciiRaster::fromFile(path);
    if (ras.ncols() != gnx || ras.nrows() != gny)
      throw std::runtime_error("RasterField: raster '" + path + "' is " +
                               std::to_string(ras.ncols()) + "x" + std::to_string(ras.nrows()) +
                               " but domain is " + std::to_string(gnx) + "x" + std::to_string(gny));
    std::vector<double> raw(static_cast<size_t>(gnx * gny));
    const RealArr1DHost& d = ras.data();
    for (int gj = 0; gj < gny; ++gj)
      for (int gi = 0; gi < gnx; ++gi) {
        const size_t src = static_cast<size_t>(gi + (gny - 1 - gj) * gnx);
        double v = d(src);
        if (ras.isNoData(static_cast<int>(src), v))
          throw std::runtime_error("RasterField: NODATA in IC raster '" + path + "'");
        raw[static_cast<size_t>(gi + gj * gnx)] = v;
      }
    return fromDoubles(raw, gnx, gny, gnz, path);
  }
  return fromDoubles(readDoublesList(path), gnx, gny, gnz, path);
}

}  // namespace frehg2
