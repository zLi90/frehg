#include "io/YamlMigration.hpp"

#include <map>
#include <sstream>
#include <string>

namespace frehg2 {

namespace {

// Move src[from] -> dst[to] only if `to` is absent and `from` is present. Removes `from`.
void renameKey(YAML::Node dst, YAML::Node src, const char* from, const char* to) {
  if (src[from] && !dst[to]) dst[to] = src[from];
  if (src[from]) src.remove(from);
}

void migrateTime(YAML::Node out) {
  // Top-level `Tend` (a common experimental spelling) lands in time.t_end.
  if (out["Tend"]) {
    if (!out["time"] || !out["time"].IsMap()) out["time"] = YAML::Node(YAML::NodeType::Map);
    if (!out["time"]["t_end"]) out["time"]["t_end"] = out["Tend"];
    out.remove("Tend");
  }
  if (!out["time"] || !out["time"].IsMap()) return;
  YAML::Node t = out["time"];
  renameKey(t, t, "Tend", "t_end");
  renameKey(t, t, "max_step", "max_steps");
  renameKey(t, t, "dt_out", "output_interval");
}

void migrateGrid(YAML::Node out) {
  if (out["grid"] && out["grid"].IsMap()) {
    if (!out["domain"] || !out["domain"].IsMap()) {
      out["domain"] = out["grid"];
    } else {
      // Merge any grid keys not already present in domain (domain wins on conflict).
      YAML::Node domain = out["domain"];
      for (const auto& kv : out["grid"]) {
        const std::string k = kv.first.as<std::string>();
        if (!domain[k]) domain[k] = kv.second;
      }
    }
    out.remove("grid");
  }
  if (out["domain"] && out["domain"].IsMap()) {
    renameKey(out["domain"], out["domain"], "bot_z", "botz");
  }
}

void migrateOutput(YAML::Node out, const std::string& sim_id) {
  if (!out["output"] || !out["output"].IsMap()) out["output"] = YAML::Node(YAML::NodeType::Map);
  YAML::Node o = out["output"];
  if (!o["filename"]) {
    std::string dir;
    if (o["directory"]) {
      dir = o["directory"].as<std::string>();
    } else if (out["io"] && out["io"]["dir"]) {
      dir = out["io"]["dir"].as<std::string>();
    }
    if (!dir.empty()) {
      if (dir.back() != '/') dir.push_back('/');
      o["filename"] = dir + sim_id + ".h5";
    }
  }
  if (o["directory"]) o.remove("directory");
  if (out["io"]) out.remove("io");
  if (!o["format"]) o["format"] = "hdf5";
}

void fixKinds(YAML::Node seq, const std::map<std::string, std::string>& renames) {
  if (!seq || !seq.IsSequence()) return;
  for (YAML::Node entry : seq) {
    if (!entry.IsMap() || !entry["type"]) continue;
    const std::string ty = entry["type"].as<std::string>();
    const auto it = renames.find(ty);
    if (it != renames.end()) entry["type"] = it->second;
  }
}

void migrateBoundaryKinds(YAML::Node out) {
  static const std::map<std::string, std::string> bc = {
      {"discharge", "bc_discharge"}, {"depth", "bc_depth"}, {"critical", "bc_critical"}};
  static const std::map<std::string, std::string> src = {{"inflow", "inflow_rate"},
                                                         {"rainfall", "rainfall_rate"},
                                                         {"well", "extraction_well"},
                                                         {"extraction", "extraction_well"}};
  if (out["boundary_conditions"] && out["boundary_conditions"].IsMap()) {
    fixKinds(out["boundary_conditions"]["surface"], bc);
    fixKinds(out["boundary_conditions"]["groundwater"], bc);
  }
  // Polygon boundaries/sources are top-level sequences in the production schema.
  fixKinds(out["boundaries"], bc);
  if (out["sources"] && out["sources"].IsSequence()) fixKinds(out["sources"], src);
}

void migrateSoil(YAML::Node out) {
  if (!out["soil"] || !out["soil"].IsMap()) return;
  YAML::Node s = out["soil"];
  if (!s["uniform"]) return;
  const bool uniform = s["uniform"].as<bool>(false);
  s.remove("uniform");
  if (!uniform || s["types"]) return;  // already has explicit classes; nothing to synthesize

  // Pull the flat uniform parameters (with frozen defaults) into a single soil class.
  auto num = [&](const char* k, double dflt) -> double {
    return s[k] ? s[k].as<double>() : dflt;
  };
  YAML::Node cls(YAML::NodeType::Map);
  cls["id"] = 0;
  cls["theta_s"] = num("theta_s", 0.4);
  cls["theta_r"] = num("theta_r", 0.0);
  YAML::Node vg(YAML::NodeType::Map);
  vg["alpha"] = num("alpha", 1.0);
  vg["n"] = num("n", 2.0);
  cls["vg"] = vg;
  YAML::Node ks(YAML::NodeType::Map);
  ks["x"] = num("ksat_x", num("ksat", 0.0));
  ks["y"] = num("ksat_y", num("ksat", 0.0));
  ks["z"] = num("ksat_z", num("ksat", 0.0));
  cls["k_sat"] = ks;

  for (const char* k : {"theta_s", "theta_r", "alpha", "n", "ksat", "ksat_x", "ksat_y", "ksat_z"}) {
    if (s[k]) s.remove(k);
  }
  YAML::Node types(YAML::NodeType::Sequence);
  types.push_back(cls);
  s["types"] = types;
  if (!s["map"]) {
    YAML::Node map(YAML::NodeType::Map);
    map["from_file"] = false;
    s["map"] = map;
  }
}

void ensureRequiredSections(YAML::Node out, const std::string& sim_id) {
  if (!out["simulation"] || !out["simulation"].IsMap()) {
    YAML::Node sim(YAML::NodeType::Map);
    sim["id"] = sim_id;
    sim["mode"] = "surface_water";
    out["simulation"] = sim;
  } else if (!out["simulation"]["id"]) {
    out["simulation"]["id"] = sim_id;
  }
  if (!out["modules"] || !out["modules"].IsMap()) {
    YAML::Node m(YAML::NodeType::Map);
    m["surface_water"] = false;
    m["groundwater"] = false;
    m["solute"] = false;
    out["modules"] = m;
  }
}

}  // namespace

YAML::Node migrateV1ToV2(const YAML::Node& in) {
  YAML::Node out = YAML::Clone(in);

  std::string sim_id = "frehg2";
  if (out["simulation"] && out["simulation"]["id"]) {
    sim_id = out["simulation"]["id"].as<std::string>();
  } else if (out["sim_id"]) {
    sim_id = out["sim_id"].as<std::string>();
  }

  migrateTime(out);
  migrateGrid(out);
  migrateOutput(out, sim_id);
  migrateBoundaryKinds(out);
  migrateSoil(out);
  ensureRequiredSections(out, sim_id);

  out["schema_version"] = "2.0";  // always normalized to the production version
  return out;
}

std::string migrateV1ToV2String(const std::string& yaml) {
  YAML::Node in = YAML::Load(yaml);
  YAML::Node out = migrateV1ToV2(in);
  std::stringstream ss;
  ss << out;
  return ss.str();
}

}  // namespace frehg2
