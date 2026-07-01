#include "bc/PolygonConfig.hpp"

#include <stdexcept>
#include <string>

#include "io/Config.hpp"

namespace frehg2 {

namespace {

Polygon parsePolygon(const YAML::Node& node, const std::string& section, size_t idx) {
  Polygon poly;
  const std::string where = section + "[" + std::to_string(idx) + "]";
  poly.name = node["name"] ? node["name"].as<std::string>() : (where);
  const YAML::Node verts = node["vertices"];
  if (!verts || !verts.IsSequence() || verts.size() < 3)
    throw std::runtime_error("PolygonConfig: " + where +
                             " ('" + poly.name + "') must have a 'vertices' list of >= 3 points");
  for (size_t k = 0; k < verts.size(); ++k) {
    const YAML::Node v = verts[k];
    if (!v.IsSequence() || v.size() != 2)
      throw std::runtime_error("PolygonConfig: " + where + " ('" + poly.name +
                               "') vertex " + std::to_string(k) + " must be [x, y]");
    try {
      poly.vertices.push_back({v[0].as<double>(), v[1].as<double>()});
    } catch (const YAML::Exception&) {
      throw std::runtime_error("PolygonConfig: " + where + " ('" + poly.name +
                               "') vertex " + std::to_string(k) + " has non-numeric coordinates");
    }
  }
  return poly;
}

std::string requireType(const YAML::Node& node, const std::string& where) {
  if (!node["type"])
    throw std::runtime_error("PolygonConfig: " + where + " is missing required key 'type'");
  return node["type"].as<std::string>();
}

real rateOf(const YAML::Node& node) {
  return node["rate"] ? node["rate"].as<double>() : 0.0;
}

}  // namespace

PolygonRegions parsePolygonRegions(const Config& cfg) {
  PolygonRegions out;
  const YAML::Node& root = cfg.root();

  const YAML::Node bnds = root["boundaries"];
  if (bnds && bnds.IsSequence()) {
    for (size_t i = 0; i < bnds.size(); ++i) {
      const YAML::Node n = bnds[i];
      const std::string where = "boundaries[" + std::to_string(i) + "]";
      const std::string type = requireType(n, where);
      BcRegion r;
      r.polygon = parsePolygon(n, "boundaries", i);
      r.kind = bcKindFromString(type);  // throws on unknown / on a source type used here
      r.rate = rateOf(n);
      out.boundaries.push_back(std::move(r));
    }
  }

  const YAML::Node srcs = root["sources"];
  // `sources:` is historically a MAP (sources.surface.rainfall...) for uniform forcing. Polygon
  // sources are a SEQUENCE. Only parse it when it is a sequence so the existing scalar-forcing
  // schema keeps working untouched.
  if (srcs && srcs.IsSequence()) {
    for (size_t i = 0; i < srcs.size(); ++i) {
      const YAML::Node n = srcs[i];
      const std::string where = "sources[" + std::to_string(i) + "]";
      const std::string type = requireType(n, where);
      SourceRegion r;
      r.polygon = parsePolygon(n, "sources", i);
      r.kind = sourceKindFromString(type);
      r.rate = rateOf(n);
      out.sources.push_back(std::move(r));
    }
  }

  return out;
}

}  // namespace frehg2
