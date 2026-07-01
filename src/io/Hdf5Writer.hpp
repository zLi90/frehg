// HDF5 implementation of OutputWriter (P3.2.1). Internal header (uses the HDF5 C API).
#ifndef FREHG2_IO_HDF5_WRITER_HPP
#define FREHG2_IO_HDF5_WRITER_HPP

#include <hdf5.h>

#include <map>
#include <set>
#include <string>
#include <vector>

#include "frehg2/io/OutputWriter.hpp"

namespace frehg2 {

class Hdf5Writer : public OutputWriter {
 public:
  Hdf5Writer(MPI_Comm comm, IoLayout layout, IoMode mode, int max_checkpoints = 2);
  ~Hdf5Writer() override;

  void openFile(const std::string& path, const RunMetadata& meta) override;
  void writeDomain(const Grid& grid) override;
  void writeSurfaceField(const std::string& name, double time,
                         const RealArr1DHost& owned_values,
                         const std::string& units) override;
  void writeSubsurfaceField(const std::string& name, double time,
                            const RealArr1DHost& owned_values,
                            const std::string& units) override;
  void writeMonitorSeries(const std::string& name, const std::vector<double>& times,
                          const std::vector<double>& values) override;
  void writeCheckpoint(const CheckpointState& cp) override;
  void close() override;

 private:
  void writeMetadata(hid_t file);
  void writeArrayWithUnits(hid_t loc, const std::string& dsname, const double* data,
                           hsize_t n, const std::string& units, bool gzip);
  // Gather owned values onto rank 0 into a global physical-ordered array.
  std::vector<double> gatherGlobal(const RealArr1DHost& owned,
                                   const IntArr1DHost& global_idx, int global_size) const;
  void writeFieldSerialGather(const std::string& group, const std::string& name,
                              double time, const RealArr1DHost& owned,
                              const IntArr1DHost& gidx, int global_size,
                              const std::string& units);
  void writeFieldParallel(const std::string& group, const std::string& name, double time,
                          const RealArr1DHost& owned, const IntArr1DHost& gidx,
                          int global_size, const std::string& units);
  void writeFieldPerRank(const std::string& group, const std::string& name, double time,
                         const RealArr1DHost& owned, const IntArr1DHost& gidx,
                         const std::string& units);
  void writeXdmf();
  std::string rankFilePath() const;

  MPI_Comm comm_;
  int rank_ = 0;
  int size_ = 1;
  IoLayout layout_;
  IoMode mode_;
  int max_checkpoints_ = 2;

  std::string path_;       // user-requested output path (the shared/canonical name)
  std::string file_stem_;  // path without extension
  RunMetadata meta_;

  hid_t file_ = -1;     // open file handle for the rank(s) that own one
  bool owns_file_ = false;
  bool grid_written_ = false;
  int grid_nx_ = 0, grid_ny_ = 0, grid_nz_ = 0;

  // For the XDMF sidecar: which (field -> {snapshot time -> dataset key}) were written.
  std::map<std::string, std::map<double, std::string>> surf_times_;
  std::map<std::string, std::map<double, std::string>> subs_times_;

  std::vector<double> checkpoint_times_;  // for pruning to max_checkpoints_
};

}  // namespace frehg2

#endif  // FREHG2_IO_HDF5_WRITER_HPP
