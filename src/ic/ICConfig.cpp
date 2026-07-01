#include "ic/ICConfig.hpp"

#include <stdexcept>

#include "io/Config.hpp"

namespace frehg2 {

namespace {

ICKind kindFromString(const std::string& s) {
  if (s == "constant") return ICKind::Constant;
  if (s == "raster" || s == "file") return ICKind::Raster;
  if (s == "polygon") return ICKind::Polygon;
  if (s == "formula") return ICKind::Formula;
  if (s == "restart") return ICKind::Restart;
  throw std::runtime_error("ICConfig: unknown initial condition type '" + s + "'");
}

ICFieldSpec parseFieldNode(const YAML::Node& node, const std::string& path, real scalar_default) {
  ICFieldSpec spec;
  spec.constant_value = scalar_default;
  // yaml-cpp: a default-constructed Node is defined+null; treat that like "missing".
  if (!node.IsDefined() || node.IsNull()) return spec;

  if (node.IsScalar()) {
    spec.kind = ICKind::Constant;
    spec.constant_value = node.as<double>();
    return spec;
  }
  if (!node.IsMap())
    throw std::runtime_error("ICConfig: '" + path + "' must be a scalar or a mapping");

  const std::string type = node["type"] ? node["type"].as<std::string>() : "constant";
  spec.kind = kindFromString(type);
  if (node["value"]) spec.constant_value = node["value"].as<double>();
  if (node["file"]) spec.file = node["file"].as<std::string>();
  if (node["format"]) spec.format = node["format"].as<std::string>();
  if (node["dataset"]) spec.dataset = node["dataset"].as<std::string>();
  if (node["formula"]) spec.formula = node["formula"].as<std::string>();
  if (node["default"]) spec.polygon_default = node["default"].as<double>();
  if (node["restart_file"]) spec.restart_file = node["restart_file"].as<std::string>();
  if (node["restart_time"]) spec.restart_time = node["restart_time"].as<double>();
  if (node["values"] && node["values"].IsSequence()) {
    for (size_t k = 0; k < node["values"].size(); ++k) {
      const YAML::Node& v = node["values"][k];
      if (!v["name"])
        throw std::runtime_error("ICConfig: polygon IC value entry missing 'name' at " + path);
      const std::string name = v["name"].as<std::string>();
      const real val = v["value"] ? v["value"].as<double>() : spec.constant_value;
      spec.polygon_values[name] = val;
    }
  }
  return spec;
}

ICFieldSpec legacyRasterIfRequested(const Config& config, const std::string& from_file_key,
                                    const std::string& file_key, ICFieldSpec base) {
  if (!config.getOr<bool>(from_file_key, false)) return base;
  ICFieldSpec spec = base;
  spec.kind = ICKind::Raster;
  spec.file = config.getOr<std::string>(file_key, "");
  if (spec.file.empty())
    throw std::runtime_error("ICConfig: " + from_file_key + " is true but " + file_key +
                             " is empty");
  return spec;
}

}  // namespace

InitialConditionsConfig parseInitialConditions(const Config& config, real default_gw_wc) {
  InitialConditionsConfig ic;
  const YAML::Node root = config.root()["initial_conditions"];
  const real default_eta = config.getOr<double>("domain.botz", 0.0);

  if (root && root["type"] && root["type"].as<std::string>() == "restart") {
    ic.use_restart = true;
    ic.restart_file = root["file"] ? root["file"].as<std::string>() : "";
    ic.restart_time = root["time"] ? root["time"].as<double>() : 0.0;
    if (ic.restart_file.empty())
      throw std::runtime_error("ICConfig: initial_conditions.type=restart requires file");
    return ic;
  }

  const YAML::Node regions = root ? root["regions"] : YAML::Node();
  if (regions && regions.IsSequence()) {
    for (size_t k = 0; k < regions.size(); ++k) {
      const YAML::Node& r = regions[k];
      if (!r["name"] || !r["vertices"])
        throw std::runtime_error("ICConfig: initial_conditions.regions entry missing name/vertices");
      ICRegion reg;
      reg.polygon.name = r["name"].as<std::string>();
      reg.value = r["value"] ? r["value"].as<double>() : 0.0;
      const YAML::Node& verts = r["vertices"];
      if (!verts.IsSequence() || verts.size() < 3)
        throw std::runtime_error("ICConfig: region '" + reg.polygon.name +
                                 "' needs at least 3 vertices");
      for (size_t vi = 0; vi < verts.size(); ++vi) {
        if (!verts[vi].IsSequence() || verts[vi].size() < 2)
          throw std::runtime_error("ICConfig: region '" + reg.polygon.name +
                                   "' vertex must be [x,y]");
        reg.polygon.vertices.push_back(
            {verts[vi][0].as<double>(), verts[vi][1].as<double>()});
      }
      ic.regions.push_back(std::move(reg));
    }
  }

  const YAML::Node sw = (root && root.IsMap()) ? root["surface_water"] : YAML::Node();
  const YAML::Node gw = (root && root.IsMap()) ? root["groundwater"] : YAML::Node();
  const YAML::Node sol = (root && root.IsMap()) ? root["solute"] : YAML::Node();

  ic.surface_eta = parseFieldNode(
      (sw && sw.IsMap()) ? sw["eta"] : YAML::Node(), "initial_conditions.surface_water.eta",
      default_eta);
  ic.surface_u = parseFieldNode((sw && sw.IsMap()) ? sw["u"] : YAML::Node(),
                                "initial_conditions.surface_water.u", 0.0);
  ic.surface_v = parseFieldNode((sw && sw.IsMap()) ? sw["v"] : YAML::Node(),
                                "initial_conditions.surface_water.v", 0.0);

  ic.groundwater_wc = parseFieldNode((gw && gw.IsMap()) ? gw["wc"] : YAML::Node(),
                                   "initial_conditions.groundwater.wc", default_gw_wc);
  ic.groundwater_head = parseFieldNode((gw && gw.IsMap()) ? gw["h"] : YAML::Node(),
                                       "initial_conditions.groundwater.h", 0.0);

  ic.solute_surface = parseFieldNode((sol && sol.IsMap()) ? sol["surface"] : YAML::Node(),
                                     "initial_conditions.solute.surface", 0.0);
  ic.solute_subsurface = parseFieldNode((sol && sol.IsMap()) ? sol["subsurface"] : YAML::Node(),
                                        "initial_conditions.solute.subsurface", 0.0);

  // Legacy *_from_file flags map to raster ICs (backward compatible with b1/b2 YAML).
  ic.surface_eta = legacyRasterIfRequested(
      config, "initial_conditions.surface_water.eta_from_file",
      "initial_conditions.surface_water.eta_file", ic.surface_eta);
  ic.groundwater_wc = legacyRasterIfRequested(
      config, "initial_conditions.groundwater.wc_from_file",
      "initial_conditions.groundwater.wc_file", ic.groundwater_wc);
  ic.groundwater_head = legacyRasterIfRequested(
      config, "initial_conditions.groundwater.h_from_file",
      "initial_conditions.groundwater.h_file", ic.groundwater_head);

  return ic;
}

}  // namespace frehg2
