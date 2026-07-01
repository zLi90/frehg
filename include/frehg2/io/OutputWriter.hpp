// Backend-agnostic output-writer interface (P3.2.0).
//
// The Orchestrator (P7) holds an `OutputWriter&` and never a concrete HDF5 type, so the
// output format is swappable (HDF5 now; netCDF/ADIOS2 later) without touching solvers.
// Fields are passed as this rank's OWNED physical cells (no halo) plus a global-index map
// (IoLayout); the writer resolves the parallel-write strategy (serial gather, parallel
// collective, or file-per-rank) so that the on-disk result always reassembles to the same
// global field.
#ifndef FREHG2_IO_OUTPUT_WRITER_HPP
#define FREHG2_IO_OUTPUT_WRITER_HPP

#include <mpi.h>

#include <memory>
#include <string>
#include <vector>

#include "frehg2/core/define.hpp"
#include "frehg2/io/Checkpoint.hpp"

namespace frehg2 {

class Grid;  // src/core/Grid.hpp (forward-declared; full type only needed in the .cpp)

// Parallel write strategy (P3.2.1a).
enum class IoMode {
  SerialGather,        // rank 0 gathers and writes one gzip-compressed file (default)
  ParallelCollective,  // all ranks write a shared file via collective MPI-IO (uncompressed)
  FilePerRank          // each rank writes its own file (+ global-index map for reassembly)
};
IoMode ioModeFromString(const std::string& s);
const char* ioModeName(IoMode m);

// Run provenance embedded in every output file (P3.2.1).
struct RunMetadata {
  std::string title = "frehg2 simulation";
  std::string version = "2.0";
  std::string date;             // ISO-8601; filled at openFile if empty
  std::string frehg2_version;
  std::string git_sha = "unknown";
  bool git_dirty = false;
  std::string config_sha256;
  std::string build_type;
  std::string compiler;
  std::string kokkos_backend;
  std::string solver_backend = "petsc";
  int mpi_ranks = 1;
};

// Maps this rank's owned cells to their global physical indices (no halo). Surface index
// is gi + gj*nx; subsurface index is gi + gj*nx + k*nx*ny (i fastest).
struct IoLayout {
  int nx = 0;
  int ny = 0;
  int nz = 0;
  IntArr1DHost surf_global_idx;  // length = owned surface cells on this rank
  IntArr1DHost subs_global_idx;  // length = owned subsurface cells on this rank
};

class OutputWriter {
 public:
  virtual ~OutputWriter() = default;

  virtual void openFile(const std::string& path, const RunMetadata& meta) = 0;
  virtual void writeDomain(const Grid& grid) = 0;

  // `owned_values` holds this rank's owned cells in the same order as the IoLayout indices.
  virtual void writeSurfaceField(const std::string& name, double time,
                                 const RealArr1DHost& owned_values,
                                 const std::string& units = "") = 0;
  virtual void writeSubsurfaceField(const std::string& name, double time,
                                    const RealArr1DHost& owned_values,
                                    const std::string& units = "") = 0;

  virtual void writeMonitorSeries(const std::string& name, const std::vector<double>& times,
                                  const std::vector<double>& values) = 0;

  virtual void writeCheckpoint(const CheckpointState& cp) = 0;

  virtual void close() = 0;
};

// Factory: selects the concrete writer by `format` ("hdf5"). The layout/comm/io_mode are
// fixed for the lifetime of the writer.
std::unique_ptr<OutputWriter> makeOutputWriter(const std::string& format, MPI_Comm comm,
                                               const IoLayout& layout, IoMode mode);

}  // namespace frehg2

#endif  // FREHG2_IO_OUTPUT_WRITER_HPP
