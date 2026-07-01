// v1/experimental -> frozen v2 YAML migration (P17, Task 17.3.2/17.3.4).
//
// The production schema is the P0-frozen "v2" schema (see docs/yaml_schema_v2.md). This module
// upgrades older/experimental YAML dialects into that frozen schema WITHOUT renaming any v2
// frozen key. The migration is purely additive/renaming on the YAML tree:
//   * time.max_step / (top-level or time.)Tend / time.dt_out  -> time.{max_steps,t_end,output_interval}
//   * grid.*  -> domain.*   (and domain.bot_z -> domain.botz)
//   * output.directory / io.dir  -> output.filename (+ output.format default)
//   * boundary/source entry `type: discharge|depth|critical|inflow|rainfall|well`
//       -> the frozen `bc_discharge|bc_depth|bc_critical|inflow_rate|rainfall_rate|extraction_well`
//   * soil.uniform: true (+ scalar soil fields) -> frozen soil.{map.from_file, types[0]} block
//   * adds schema_version: "2.0" and any missing required top-level sections with sane defaults
//
// Migrating an already-v2 YAML is idempotent (the renames only fire when a v1 key is present and
// the corresponding frozen key is absent), so all benchmark v2 fixtures round-trip unchanged.
#ifndef FREHG2_IO_YAML_MIGRATION_HPP
#define FREHG2_IO_YAML_MIGRATION_HPP

#include <yaml-cpp/yaml.h>

#include <string>

namespace frehg2 {

// Migrate a parsed YAML tree to the frozen v2 schema. Returns a new tree (input is not mutated).
YAML::Node migrateV1ToV2(const YAML::Node& in);

// Convenience: parse `yaml`, migrate, and emit the v2 YAML as a string.
std::string migrateV1ToV2String(const std::string& yaml);

}  // namespace frehg2

#endif  // FREHG2_IO_YAML_MIGRATION_HPP
