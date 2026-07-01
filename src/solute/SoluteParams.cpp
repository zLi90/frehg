#include "solute/SoluteParams.hpp"

#include "io/Config.hpp"

namespace frehg2 {

SoluteParams SoluteParams::fromConfig(const Config& cfg) {
  SoluteParams p;
  p.enabled = cfg.getOr<bool>("solute.enabled", p.enabled);
  p.c_rain = cfg.getOr<double>("solute.c_rain", p.c_rain);
  p.k_decay = cfg.getOr<double>("solute.k_decay", p.k_decay);
  p.D = cfg.getOr<double>("solute.D", p.D);
  p.advection_scheme = cfg.getOr<std::string>("solute.advection_scheme", p.advection_scheme);
  p.diffusion_scheme = cfg.getOr<std::string>("solute.diffusion_scheme", p.diffusion_scheme);
  p.cfl_max = cfg.getOr<double>("solute.cfl_max", p.cfl_max);
  p.min_depth = cfg.getOr<double>("solute.min_depth", p.min_depth);
  return p;
}

}  // namespace frehg2
