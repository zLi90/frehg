// HDF5 reader (P3.2.2 / P3.5.2). Reads fields and checkpoints back as host arrays.
#ifndef FREHG2_IO_HDF5_READER_HPP
#define FREHG2_IO_HDF5_READER_HPP

#include <string>
#include <vector>

#include "frehg2/io/Checkpoint.hpp"
#include "frehg2/core/define.hpp"

namespace frehg2 {

class Hdf5Reader {
 public:
  // Read a field written in serial_gather / parallel_collective mode (single shared file).
  // `time_step` is the integer time key used at write time.
  static RealArr1DHost readSurface(const std::string& path, const std::string& field,
                                   long long time_step);
  static RealArr1DHost readSubsurface(const std::string& path, const std::string& field,
                                      long long time_step);

  // Reassemble a field written in file_per_rank mode from all rank files, using the stored
  // global-index maps. `nranks` rank files of the form <stem>.rankN.h5 are read.
  static RealArr1DHost readSurfaceFilePerRank(const std::string& stem,
                                              const std::string& field, long long time_step,
                                              int nranks, int global_size);

  // Read a string attribute from /simulation/metadata.
  static std::string readMetadataAttr(const std::string& path, const std::string& attr);

  // Read the per-field `units` attribute on a surface dataset.
  static std::string readSurfaceUnits(const std::string& path, const std::string& field,
                                      long long time_step);

  // Read a checkpoint file. If expected_config_sha256 is non-empty, the returned
  // RestartState.config_matches reflects whether it matched.
  static RestartState readCheckpoint(const std::string& path,
                                     const std::string& expected_config_sha256 = "");
};

}  // namespace frehg2

#endif  // FREHG2_IO_HDF5_READER_HPP
