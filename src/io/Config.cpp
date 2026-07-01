#include "io/Config.hpp"

#include <sstream>

namespace frehg2 {

namespace {
std::vector<std::string> splitDot(const std::string& path) {
  std::vector<std::string> keys;
  std::string cur;
  for (char c : path) {
    if (c == '.') {
      keys.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  keys.push_back(cur);
  return keys;
}

std::string dirOf(const std::string& path) {
  const auto slash = path.find_last_of('/');
  if (slash == std::string::npos) return "";
  return path.substr(0, slash);
}

// Recursive const descent. Using the *const* Node::operator[] is important: it returns an
// undefined node for a missing key WITHOUT auto-creating it, and never mutates the tree
// (unlike `node = node[key]` on a non-const alias, which assigns content into the alias).
YAML::Node descend(const YAML::Node& node, const std::vector<std::string>& keys, size_t i,
                   bool& found) {
  if (i == keys.size()) {
    found = true;
    return node;
  }
  if (!node.IsMap()) {
    found = false;
    return YAML::Node();
  }
  const YAML::Node child = node[keys[i]];
  if (!child.IsDefined() || child.IsNull()) {
    found = false;
    return YAML::Node();
  }
  return descend(child, keys, i + 1, found);
}
}  // namespace

Config Config::fromFile(const std::string& path) {
  Config cfg;
  cfg.source_ = path;
  cfg.config_dir_ = dirOf(path);
  try {
    cfg.root_ = YAML::LoadFile(path);
  } catch (const YAML::Exception& e) {
    throw std::runtime_error("failed to parse YAML file '" + path + "': " + e.what());
  }
  cfg.validateSchema();
  return cfg;
}

Config Config::fromString(const std::string& yaml, const std::string& config_dir) {
  Config cfg;
  cfg.source_ = "<string>";
  cfg.config_dir_ = config_dir;
  try {
    cfg.root_ = YAML::Load(yaml);
  } catch (const YAML::Exception& e) {
    throw std::runtime_error(std::string("failed to parse YAML string: ") + e.what());
  }
  cfg.validateSchema();
  return cfg;
}

YAML::Node Config::findNode(const std::string& path, bool& found) const {
  found = false;
  return descend(root_, splitDot(path), 0, found);
}

size_t Config::indexedCount(const std::string& path) const {
  bool found = false;
  YAML::Node n = findNode(path, found);
  if (!found || !n.IsSequence()) return 0;
  return n.size();
}

std::string Config::resolvePath(const std::string& rel) const {
  if (rel.empty() || rel[0] == '/') return rel;            // absolute or empty
  if (config_dir_.empty()) return rel;                      // relative to cwd
  return config_dir_ + "/" + rel;
}

std::string Config::resolvedString() const {
  std::stringstream ss;
  ss << root_;
  return ss.str();
}

void Config::validateSchema() const {
  // modules section present
  bool found = false;
  findNode("modules", found);
  if (!found) {
    throw std::runtime_error("schema error: 'modules' section missing (in " + source_ + ")");
  }
  // time.t_end and time.dt present and numeric
  for (const char* key : {"time.t_end", "time.dt"}) {
    bool f = false;
    YAML::Node n = findNode(key, f);
    if (!f) {
      throw std::runtime_error(std::string("schema error: required key '") + key +
                               "' missing (in " + source_ + ")");
    }
    try {
      (void)n.as<double>();
    } catch (const YAML::Exception&) {
      throw std::runtime_error(std::string("schema error: '") + key +
                               "' must be numeric (in " + source_ + ")");
    }
  }
  // domain.nx/ny/nz present and positive integers
  for (const char* key : {"domain.nx", "domain.ny", "domain.nz"}) {
    bool f = false;
    YAML::Node n = findNode(key, f);
    if (!f) {
      throw std::runtime_error(std::string("schema error: required key '") + key +
                               "' missing (in " + source_ + ")");
    }
    int v = 0;
    try {
      v = n.as<int>();
    } catch (const YAML::Exception&) {
      throw std::runtime_error(std::string("schema error: '") + key +
                               "' must be an integer (in " + source_ + ")");
    }
    if (v <= 0) {
      throw std::runtime_error(std::string("schema error: '") + key +
                               "' must be positive (in " + source_ + ")");
    }
  }
}

}  // namespace frehg2
