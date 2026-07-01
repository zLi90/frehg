#include "monitoring/MonitorConfig.hpp"

#include <stdexcept>

#include "io/Config.hpp"

namespace frehg2 {

namespace {

std::vector<std::string> defaultSwFields() { return {"eta", "depth", "u", "v"}; }
std::vector<std::string> defaultGwFields() { return {"head", "moisture"}; }

std::vector<std::string> parseFieldList(const YAML::Node& node) {
  std::vector<std::string> out;
  if (!node || !node.IsSequence()) return out;
  for (size_t k = 0; k < node.size(); ++k) out.push_back(node[k].as<std::string>());
  return out;
}

void parseProbeNode(const YAML::Node& node, const std::string& where, bool subsurface,
                    std::vector<std::string> default_fields, MonitorBundle& out) {
  ProbeSpec p;
  p.subsurface = subsurface;
  p.name = node["name"] ? node["name"].as<std::string>() : where;
  if (node["xyz"] && node["xyz"].IsSequence() && node["xyz"].size() >= 2) {
    p.coords_from_xyz = true;
    p.x = node["xyz"][0].as<double>();
    p.y = node["xyz"][1].as<double>();
    if (node["xyz"].size() >= 3) p.z = node["xyz"][2].as<double>();
  } else {
    if (!node["i"] || !node["j"])
      throw std::runtime_error("MonitorConfig: probe '" + p.name +
                               "' needs xyz or i/j indices");
    p.indices_from_grid = true;
    p.gi = node["i"].as<int>();
    p.gj = node["j"].as<int>();
    if (node["k"]) p.gk = node["k"].as<int>();
    if (subsurface && !node["k"]) p.gk = 0;
  }
  p.fields = parseFieldList(node["fields"]);
  if (p.fields.empty()) p.fields = std::move(default_fields);
  out.probes.push_back(std::move(p));
}

void parseLineNode(const YAML::Node& node, const std::string& where, bool subsurface,
                   MonitorBundle& out) {
  LineFluxSpec line;
  line.subsurface = subsurface;
  line.name = node["name"] ? node["name"].as<std::string>() : where;
  if (node["p0"] && node["p0"].IsSequence() && node["p0"].size() >= 2) {
    line.p0[0] = node["p0"][0].as<double>();
    line.p0[1] = node["p0"][1].as<double>();
    if (node["p0"].size() >= 3) line.p0[2] = node["p0"][2].as<double>();
  }
  if (node["p1"] && node["p1"].IsSequence() && node["p1"].size() >= 2) {
    line.p1[0] = node["p1"][0].as<double>();
    line.p1[1] = node["p1"][1].as<double>();
    if (node["p1"].size() >= 3) line.p1[2] = node["p1"][2].as<double>();
  } else {
    throw std::runtime_error("MonitorConfig: line '" + line.name + "' requires p0 and p1");
  }
  if (node["field"]) line.field = node["field"].as<std::string>();
  out.lines.push_back(std::move(line));
}

void parseProbesSequence(const YAML::Node& seq, bool subsurface,
                         const std::vector<std::string>& default_fields, MonitorBundle& out) {
  if (!seq || !seq.IsSequence()) return;
  for (size_t k = 0; k < seq.size(); ++k) {
    const std::string where =
        std::string(subsurface ? "monitors.probes_subsurface" : "monitors.probes") + "[" +
        std::to_string(k) + "]";
    parseProbeNode(seq[k], where, subsurface, default_fields, out);
  }
}

void parseLinesSequence(const YAML::Node& seq, bool subsurface, MonitorBundle& out) {
  if (!seq || !seq.IsSequence()) return;
  for (size_t k = 0; k < seq.size(); ++k) {
    const std::string where =
        std::string(subsurface ? "monitors.lines_subsurface" : "monitors.lines") + "[" +
        std::to_string(k) + "]";
    parseLineNode(seq[k], where, subsurface, out);
  }
}

}  // namespace

MonitorBundle parseMonitors(const Config& config, bool sw_enabled, bool gw_enabled) {
  MonitorBundle bundle;
  const YAML::Node root = config.root();

  const YAML::Node monitors = root["monitors"];
  if (monitors && monitors.IsMap()) {
    parseProbesSequence(monitors["probes"], false, sw_enabled ? defaultSwFields() : std::vector<std::string>{},
                        bundle);
    parseProbesSequence(monitors["probes_subsurface"], true,
                        gw_enabled ? defaultGwFields() : std::vector<std::string>{}, bundle);
    parseLinesSequence(monitors["lines"], false, bundle);
    parseLinesSequence(monitors["lines_subsurface"], true, bundle);
  }

  // Legacy / benchmark schema-2.0: monitoring.points (grid indices).
  const YAML::Node legacy = root["monitoring"];
  if (legacy && legacy.IsMap()) {
    const YAML::Node points = legacy["points"];
    if (points && points.IsSequence()) {
      for (size_t k = 0; k < points.size(); ++k) {
        const YAML::Node& pt = points[k];
        ProbeSpec p;
        p.name = pt["name"] ? pt["name"].as<std::string>()
                              : ("monitor_" + std::to_string(k));
        if (!pt["i"] || !pt["j"])
          throw std::runtime_error("MonitorConfig: monitoring.points[" + std::to_string(k) +
                                   "] requires i and j");
        p.gi = pt["i"].as<int>();
        p.gj = pt["j"].as<int>();
        p.indices_from_grid = true;
        if (pt["k"]) {
          p.gk = pt["k"].as<int>();
          p.subsurface = true;
        } else {
          p.subsurface = gw_enabled && !sw_enabled;
        }
        p.fields = parseFieldList(pt["fields"]);
        if (p.fields.empty()) {
          p.fields = p.subsurface ? defaultGwFields() : defaultSwFields();
        }
        bundle.probes.push_back(std::move(p));
      }
    }
  }

  return bundle;
}

}  // namespace frehg2
