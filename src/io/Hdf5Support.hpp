// Internal HDF5 C-API helpers + RAII handle guard (P3.2). Shared by Hdf5Writer/Reader.
// Not part of the public interface. Uses the HDF5 **C** API (no H5Cpp.h locally).
#ifndef FREHG2_IO_HDF5_SUPPORT_HPP
#define FREHG2_IO_HDF5_SUPPORT_HPP

#include <hdf5.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace frehg2 {
namespace h5 {

// RAII for an HDF5 object id. All H5?close functions share the herr_t(hid_t) signature.
class Guard {
 public:
  Guard() = default;
  Guard(hid_t id, herr_t (*closer)(hid_t)) : id_(id), closer_(closer) {}
  ~Guard() { reset(); }
  Guard(const Guard&) = delete;
  Guard& operator=(const Guard&) = delete;
  Guard(Guard&& o) noexcept : id_(o.id_), closer_(o.closer_) {
    o.id_ = -1;
    o.closer_ = nullptr;
  }
  Guard& operator=(Guard&& o) noexcept {
    if (this != &o) {
      reset();
      id_ = o.id_;
      closer_ = o.closer_;
      o.id_ = -1;
      o.closer_ = nullptr;
    }
    return *this;
  }
  void reset() {
    if (id_ >= 0 && closer_) closer_(id_);
    id_ = -1;
    closer_ = nullptr;
  }
  hid_t get() const { return id_; }
  bool valid() const { return id_ >= 0; }

 private:
  hid_t id_ = -1;
  herr_t (*closer_)(hid_t) = nullptr;
};

inline void check(hid_t id, const char* what) {
  if (id < 0) throw std::runtime_error(std::string("HDF5 error: ") + what);
}
inline void checkErr(herr_t e, const char* what) {
  if (e < 0) throw std::runtime_error(std::string("HDF5 error: ") + what);
}

// Create all groups along a '/'-separated path (idempotent); returns the leaf group.
Guard ensureGroup(hid_t file, const std::string& path);

bool linkExists(hid_t loc, const std::string& path);

// Write a 1D double dataset at loc/name (optionally gzip-compressed + chunked).
void writeDoubleDataset(hid_t loc, const std::string& name, const double* data, hsize_t n,
                        bool gzip);
void writeIntDataset(hid_t loc, const std::string& name, const int* data, hsize_t n);

// Read a 1D dataset (any numeric type) at loc/name into a double vector.
std::vector<double> readDoubleDataset(hid_t loc, const std::string& name);
std::vector<int> readIntDataset(hid_t loc, const std::string& name);

void writeStringAttr(hid_t loc, const std::string& name, const std::string& value);
void writeDoubleAttr(hid_t loc, const std::string& name, double value);
void writeIntAttr(hid_t loc, const std::string& name, int value);

std::string readStringAttr(hid_t loc, const std::string& name);
double readDoubleAttr(hid_t loc, const std::string& name);
int readIntAttr(hid_t loc, const std::string& name);
bool attrExists(hid_t loc, const std::string& name);

}  // namespace h5
}  // namespace frehg2

#endif  // FREHG2_IO_HDF5_SUPPORT_HPP
