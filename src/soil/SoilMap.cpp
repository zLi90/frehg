#include "soil/SoilMap.hpp"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace frehg2 {

void SoilMap::setClasses(std::vector<SoilParams> classes) {
  if (classes.empty())
    throw std::runtime_error("SoilMap::setClasses: class table must have at least one class");
  classes_ = std::move(classes);
}

void SoilMap::setClassIndex(int nx, int ny, std::vector<int> class_idx) {
  if (classes_.empty())
    throw std::runtime_error("SoilMap::setClassIndex: set classes before the class index");
  if (nx <= 0 || ny <= 0)
    throw std::runtime_error("SoilMap::setClassIndex: nx and ny must be positive");
  const size_t expect = static_cast<size_t>(nx) * static_cast<size_t>(ny);
  if (class_idx.size() != expect)
    throw std::runtime_error("SoilMap::setClassIndex: got " + std::to_string(class_idx.size()) +
                             " indices, expected nx*ny=" + std::to_string(expect));
  const int nc = nClasses();
  for (size_t k = 0; k < class_idx.size(); ++k) {
    if (class_idx[k] < 0 || class_idx[k] >= nc)
      throw std::runtime_error("SoilMap::setClassIndex: class id " + std::to_string(class_idx[k]) +
                               " at cell " + std::to_string(k) + " is out of range [0," +
                               std::to_string(nc) + ")");
  }
  nx_ = nx;
  ny_ = ny;
  nz_ = 0;
  class_idx_ = std::move(class_idx);
}

void SoilMap::setClassIndex3D(int nx, int ny, int nz, std::vector<int> class_idx) {
  if (classes_.empty())
    throw std::runtime_error("SoilMap::setClassIndex3D: set classes before the class index");
  if (nx <= 0 || ny <= 0 || nz <= 0)
    throw std::runtime_error("SoilMap::setClassIndex3D: nx, ny and nz must be positive");
  const size_t expect =
      static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz);
  if (class_idx.size() != expect)
    throw std::runtime_error("SoilMap::setClassIndex3D: got " + std::to_string(class_idx.size()) +
                             " indices, expected nx*ny*nz=" + std::to_string(expect));
  const int nc = nClasses();
  for (size_t c = 0; c < class_idx.size(); ++c) {
    if (class_idx[c] < 0 || class_idx[c] >= nc)
      throw std::runtime_error("SoilMap::setClassIndex3D: class id " + std::to_string(class_idx[c]) +
                               " at cell " + std::to_string(c) + " is out of range [0," +
                               std::to_string(nc) + ")");
  }
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  class_idx_ = std::move(class_idx);
}

void SoilMap::setUniform(int nx, int ny, const SoilParams& cls) {
  setClasses({cls});
  setClassIndex(nx, ny, std::vector<int>(static_cast<size_t>(nx) * static_cast<size_t>(ny), 0));
}

void SoilMap::loadIndexFromCSV(const std::string& path, int nx, int ny) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("SoilMap::loadIndexFromCSV: cannot open '" + path + "'");
  const size_t expect = static_cast<size_t>(nx) * static_cast<size_t>(ny);
  std::vector<int> idx;
  idx.reserve(expect);
  // Treat commas as whitespace so plain CSV and whitespace-separated files both parse.
  std::string token;
  auto flush = [&](std::string& t) {
    if (t.empty()) return;
    try {
      idx.push_back(std::stoi(t));
    } catch (const std::exception&) {
      throw std::runtime_error("SoilMap::loadIndexFromCSV: non-integer token '" + t + "' in '" +
                               path + "'");
    }
    t.clear();
  };
  char ch;
  while (in.get(ch)) {
    if (ch == ',' || std::isspace(static_cast<unsigned char>(ch))) {
      flush(token);
    } else {
      token.push_back(ch);
    }
  }
  flush(token);
  if (idx.size() != expect)
    throw std::runtime_error("SoilMap::loadIndexFromCSV: '" + path + "' has " +
                             std::to_string(idx.size()) + " values, expected nx*ny=" +
                             std::to_string(expect));
  setClassIndex(nx, ny, std::move(idx));
}

bool SoilMap::isUniform() const {
  if (empty()) return true;
  const int id0 = class_idx_[0];
  for (int v : class_idx_)
    if (v != id0) return false;
  return true;
}

}  // namespace frehg2
