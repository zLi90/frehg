// Checkpoint / restart state containers (P3.5). Backend-neutral (no HDF5 here).
#ifndef FREHG2_IO_CHECKPOINT_HPP
#define FREHG2_IO_CHECKPOINT_HPP

#include <map>
#include <string>

#include "frehg2/core/define.hpp"

namespace frehg2 {

// Solver state captured at a checkpoint. Arrays are stored as given (the layout string
// records whether they are halo-backed internal storage or physical-no-halo).
struct CheckpointState {
  double time = 0.0;
  double dt = 0.0;

  bool has_sw = false;
  bool has_gw = false;
  bool has_solute = false;

  // Surface water (halo-backed internal storage by default).
  RealArr1DHost eta;
  RealArr1DHost u;
  RealArr1DHost v;

  // Groundwater.
  RealArr1DHost h;
  RealArr1DHost hn;
  RealArr1DHost wc;
  RealArr1DHost wcn;

  // Solute (if active).
  RealArr1DHost C;

  // Complete, backend-neutral solver state for bit-exact restart (P7.4). The Orchestrator
  // dumps EVERY halo-padded solver field here (keys prefixed "sw_"/"gw_") so restart needs
  // no field re-derivation and avoids the solvers' one-step geometry/drag lag. The named
  // eta/u/v/h/... members above stay populated for the documented schema + h5py interop.
  std::map<std::string, RealArr1DHost> extra;

  std::string storage_layout = "internal-with-halo";  // or "physical-no-halo"
  std::string config_sha256;
  std::string git_sha;
};

// Result of reading a checkpoint back.
struct RestartState {
  double time = 0.0;
  double dt = 0.0;

  bool has_sw = false;
  bool has_gw = false;
  bool has_solute = false;

  RealArr1DHost eta;
  RealArr1DHost u;
  RealArr1DHost v;
  RealArr1DHost h;
  RealArr1DHost hn;
  RealArr1DHost wc;
  RealArr1DHost wcn;
  RealArr1DHost C;

  // Full solver state restored from "checkpoint/extra" (see CheckpointState::extra).
  std::map<std::string, RealArr1DHost> extra;

  std::string storage_layout;
  std::string config_sha256;
  std::string git_sha;

  // True if config_sha256 matched the value the caller expected (set by readCheckpoint).
  bool config_matches = true;
};

}  // namespace frehg2

#endif  // FREHG2_IO_CHECKPOINT_HPP
