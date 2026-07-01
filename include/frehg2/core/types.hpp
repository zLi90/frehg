// Frehg2 domain/time/module parameter structures (P1.3).
//
// Placed under include/frehg2/core/ per the .cursorrules "Header files: include/frehg2/**"
// rule (the plan text's src/core/types.hpp is superseded by the public-header layout).
#ifndef FREHG2_CORE_TYPES_HPP
#define FREHG2_CORE_TYPES_HPP

#include "frehg2/core/define.hpp"

namespace frehg2 {

// Grid geometry. `nz`, `dz`, `dz_incre`, `bot_z` are subsurface (GW) parameters; the
// vertical layering is geometric (dz[k] = dz[k-1] * dz_incre), k=0 = top active layer
// (see docs/legacy_audit/index_conventions.md).
struct DomainParams {
  int nx = 1;
  int ny = 1;
  int nz = 1;
  real dx = 1.0;
  real dy = 1.0;
  real dz = 1.0;
  real dz_incre = 1.0;
  real bot_z = 0.0;
};

// Time-stepping and output/checkpoint cadence.
struct TimeParams {
  real dt = 0.0;
  real t_end = 0.0;
  real output_interval = 0.0;
  int max_steps = 0;
  real dt_checkpoint = 0.0;
  int max_checkpoints = 0;
};

// Which physics modules are enabled for a run.
struct ModuleFlags {
  bool surface_water = false;
  bool groundwater = false;
  bool solute = false;
};

}  // namespace frehg2

#endif  // FREHG2_CORE_TYPES_HPP
