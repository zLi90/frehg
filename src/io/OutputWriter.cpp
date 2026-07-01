#include "frehg2/io/OutputWriter.hpp"

#include <stdexcept>

#include "io/Hdf5Writer.hpp"

namespace frehg2 {

IoMode ioModeFromString(const std::string& s) {
  if (s == "serial_gather") return IoMode::SerialGather;
  if (s == "parallel_collective") return IoMode::ParallelCollective;
  if (s == "file_per_rank") return IoMode::FilePerRank;
  throw std::runtime_error("unknown output.io_mode: '" + s +
                           "' (expected serial_gather|parallel_collective|file_per_rank)");
}

const char* ioModeName(IoMode m) {
  switch (m) {
    case IoMode::SerialGather:
      return "serial_gather";
    case IoMode::ParallelCollective:
      return "parallel_collective";
    case IoMode::FilePerRank:
      return "file_per_rank";
  }
  return "serial_gather";
}

std::unique_ptr<OutputWriter> makeOutputWriter(const std::string& format, MPI_Comm comm,
                                               const IoLayout& layout, IoMode mode) {
  if (format == "hdf5") {
    return std::make_unique<Hdf5Writer>(comm, layout, mode);
  }
  throw std::runtime_error("unknown output.format: '" + format + "' (only 'hdf5' for now)");
}

}  // namespace frehg2
