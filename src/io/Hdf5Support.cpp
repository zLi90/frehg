#include "io/Hdf5Support.hpp"

#include <cstring>

namespace frehg2 {
namespace h5 {

namespace {
std::vector<std::string> splitSlash(const std::string& path) {
  std::vector<std::string> parts;
  std::string cur;
  for (char c : path) {
    if (c == '/') {
      if (!cur.empty()) parts.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  if (!cur.empty()) parts.push_back(cur);
  return parts;
}
}  // namespace

Guard ensureGroup(hid_t file, const std::string& path) {
  hid_t loc = file;
  Guard last;
  for (const auto& part : splitSlash(path)) {
    if (H5Lexists(loc, part.c_str(), H5P_DEFAULT) > 0) {
      hid_t g = H5Gopen2(loc, part.c_str(), H5P_DEFAULT);
      check(g, ("H5Gopen " + part).c_str());
      last = Guard(g, H5Gclose);
    } else {
      hid_t g = H5Gcreate2(loc, part.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      check(g, ("H5Gcreate " + part).c_str());
      last = Guard(g, H5Gclose);
    }
    loc = last.get();
  }
  return last;
}

bool linkExists(hid_t loc, const std::string& path) {
  // Walk the path so intermediate-missing links don't trigger HDF5 errors.
  hid_t cur = loc;
  Guard holder;
  for (const auto& part : splitSlash(path)) {
    if (H5Lexists(cur, part.c_str(), H5P_DEFAULT) <= 0) return false;
    H5O_info2_t info;
    if (H5Oget_info_by_name3(cur, part.c_str(), &info, H5O_INFO_BASIC, H5P_DEFAULT) < 0) {
      return false;
    }
    if (info.type == H5O_TYPE_GROUP) {
      hid_t g = H5Gopen2(cur, part.c_str(), H5P_DEFAULT);
      if (g < 0) return false;
      holder = Guard(g, H5Gclose);
      cur = g;
    }
  }
  return true;
}

void writeDoubleDataset(hid_t loc, const std::string& name, const double* data, hsize_t n,
                        bool gzip) {
  hsize_t dims[1] = {n};
  Guard space(H5Screate_simple(1, dims, nullptr), H5Sclose);
  check(space.get(), "H5Screate_simple");

  Guard dcpl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose);
  check(dcpl.get(), "H5Pcreate dcpl");
  if (gzip && n > 0) {
    hsize_t chunk[1] = {n};
    checkErr(H5Pset_chunk(dcpl.get(), 1, chunk), "H5Pset_chunk");
    checkErr(H5Pset_deflate(dcpl.get(), 6), "H5Pset_deflate");
  }
  Guard dset(H5Dcreate2(loc, name.c_str(), H5T_IEEE_F64LE, space.get(), H5P_DEFAULT,
                        dcpl.get(), H5P_DEFAULT),
             H5Dclose);
  check(dset.get(), ("H5Dcreate " + name).c_str());
  if (n > 0) {
    checkErr(H5Dwrite(dset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data),
             ("H5Dwrite " + name).c_str());
  }
}

void writeIntDataset(hid_t loc, const std::string& name, const int* data, hsize_t n) {
  hsize_t dims[1] = {n};
  Guard space(H5Screate_simple(1, dims, nullptr), H5Sclose);
  check(space.get(), "H5Screate_simple int");
  Guard dset(H5Dcreate2(loc, name.c_str(), H5T_STD_I32LE, space.get(), H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT),
             H5Dclose);
  check(dset.get(), ("H5Dcreate int " + name).c_str());
  if (n > 0) {
    checkErr(H5Dwrite(dset.get(), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data),
             ("H5Dwrite int " + name).c_str());
  }
}

std::vector<double> readDoubleDataset(hid_t loc, const std::string& name) {
  Guard dset(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), H5Dclose);
  check(dset.get(), ("H5Dopen " + name).c_str());
  Guard space(H5Dget_space(dset.get()), H5Sclose);
  const hssize_t n = H5Sget_simple_extent_npoints(space.get());
  std::vector<double> out(static_cast<size_t>(n));
  if (n > 0) {
    checkErr(H5Dread(dset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     out.data()),
             ("H5Dread " + name).c_str());
  }
  return out;
}

std::vector<int> readIntDataset(hid_t loc, const std::string& name) {
  Guard dset(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), H5Dclose);
  check(dset.get(), ("H5Dopen int " + name).c_str());
  Guard space(H5Dget_space(dset.get()), H5Sclose);
  const hssize_t n = H5Sget_simple_extent_npoints(space.get());
  std::vector<int> out(static_cast<size_t>(n));
  if (n > 0) {
    checkErr(H5Dread(dset.get(), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()),
             ("H5Dread int " + name).c_str());
  }
  return out;
}

void writeStringAttr(hid_t loc, const std::string& name, const std::string& value) {
  Guard atype(H5Tcopy(H5T_C_S1), H5Tclose);
  checkErr(H5Tset_size(atype.get(), value.size() + 1), "H5Tset_size");
  checkErr(H5Tset_strpad(atype.get(), H5T_STR_NULLTERM), "H5Tset_strpad");
  Guard space(H5Screate(H5S_SCALAR), H5Sclose);
  if (H5Aexists(loc, name.c_str()) > 0) H5Adelete(loc, name.c_str());
  Guard attr(H5Acreate2(loc, name.c_str(), atype.get(), space.get(), H5P_DEFAULT,
                        H5P_DEFAULT),
             H5Aclose);
  check(attr.get(), ("H5Acreate str " + name).c_str());
  checkErr(H5Awrite(attr.get(), atype.get(), value.c_str()), ("H5Awrite str " + name).c_str());
}

void writeDoubleAttr(hid_t loc, const std::string& name, double value) {
  Guard space(H5Screate(H5S_SCALAR), H5Sclose);
  if (H5Aexists(loc, name.c_str()) > 0) H5Adelete(loc, name.c_str());
  Guard attr(H5Acreate2(loc, name.c_str(), H5T_IEEE_F64LE, space.get(), H5P_DEFAULT,
                        H5P_DEFAULT),
             H5Aclose);
  check(attr.get(), ("H5Acreate dbl " + name).c_str());
  checkErr(H5Awrite(attr.get(), H5T_NATIVE_DOUBLE, &value), ("H5Awrite dbl " + name).c_str());
}

void writeIntAttr(hid_t loc, const std::string& name, int value) {
  Guard space(H5Screate(H5S_SCALAR), H5Sclose);
  if (H5Aexists(loc, name.c_str()) > 0) H5Adelete(loc, name.c_str());
  Guard attr(H5Acreate2(loc, name.c_str(), H5T_STD_I32LE, space.get(), H5P_DEFAULT,
                        H5P_DEFAULT),
             H5Aclose);
  check(attr.get(), ("H5Acreate int " + name).c_str());
  checkErr(H5Awrite(attr.get(), H5T_NATIVE_INT, &value), ("H5Awrite int " + name).c_str());
}

bool attrExists(hid_t loc, const std::string& name) {
  return H5Aexists(loc, name.c_str()) > 0;
}

std::string readStringAttr(hid_t loc, const std::string& name) {
  Guard attr(H5Aopen(loc, name.c_str(), H5P_DEFAULT), H5Aclose);
  check(attr.get(), ("H5Aopen str " + name).c_str());
  Guard atype(H5Aget_type(attr.get()), H5Tclose);
  const size_t sz = H5Tget_size(atype.get());
  std::vector<char> buf(sz + 1, '\0');
  checkErr(H5Aread(attr.get(), atype.get(), buf.data()), ("H5Aread str " + name).c_str());
  return std::string(buf.data());
}

double readDoubleAttr(hid_t loc, const std::string& name) {
  Guard attr(H5Aopen(loc, name.c_str(), H5P_DEFAULT), H5Aclose);
  check(attr.get(), ("H5Aopen dbl " + name).c_str());
  double v = 0.0;
  checkErr(H5Aread(attr.get(), H5T_NATIVE_DOUBLE, &v), ("H5Aread dbl " + name).c_str());
  return v;
}

int readIntAttr(hid_t loc, const std::string& name) {
  Guard attr(H5Aopen(loc, name.c_str(), H5P_DEFAULT), H5Aclose);
  check(attr.get(), ("H5Aopen int " + name).c_str());
  int v = 0;
  checkErr(H5Aread(attr.get(), H5T_NATIVE_INT, &v), ("H5Aread int " + name).c_str());
  return v;
}

}  // namespace h5
}  // namespace frehg2
