// YAML configuration reader (P3.1).
//
// Thin, fail-loud wrapper over yaml-cpp with dotted-path access ("time.t_end"), relative
// file-path resolution against the config's directory, and minimal schema validation.
// Missing required keys throw with the offending key and file name; nothing is silently
// defaulted (use getOr for optional keys).
#ifndef FREHG2_IO_CONFIG_HPP
#define FREHG2_IO_CONFIG_HPP

#include <yaml-cpp/yaml.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace frehg2 {

class Config {
 public:
  // Load from a YAML file; resolves relative paths against the file's directory and runs
  // schema validation. Throws std::runtime_error on parse/validation failure.
  static Config fromFile(const std::string& path);
  // Load from an in-memory YAML string (config_dir used for path resolution; "" = cwd).
  static Config fromString(const std::string& yaml, const std::string& config_dir = "");

  // Required key: throws "missing required key: <path> (in <file>)" if absent, or
  // "wrong type for key: <path>" if the value cannot convert to T.
  template <typename T>
  T get(const std::string& path) const {
    bool found = false;
    YAML::Node n = findNode(path, found);
    if (!found) {
      throw std::runtime_error("missing required key: " + path + " (in " + source_ + ")");
    }
    try {
      return n.as<T>();
    } catch (const YAML::Exception&) {
      throw std::runtime_error("wrong type for key: " + path + " (in " + source_ + ")");
    }
  }

  // Optional key: returns default_val if absent; throws on wrong type if present.
  template <typename T>
  T getOr(const std::string& path, const T& default_val) const {
    bool found = false;
    YAML::Node n = findNode(path, found);
    if (!found) return default_val;
    try {
      return n.as<T>();
    } catch (const YAML::Exception&) {
      throw std::runtime_error("wrong type for key: " + path + " (in " + source_ + ")");
    }
  }

  bool has(const std::string& path) const {
    bool found = false;
    findNode(path, found);
    return found;
  }

  // Number of items in a YAML sequence at path (0 if missing or not a sequence).
  size_t indexedCount(const std::string& path) const;

  // Resolve a (possibly relative) path against the config directory.
  std::string resolvePath(const std::string& rel) const;

  const std::string& source() const { return source_; }
  const std::string& configDir() const { return config_dir_; }

  // Canonical serialization of the resolved config (used for the provenance hash).
  std::string resolvedString() const;

  const YAML::Node& root() const { return root_; }

 private:
  Config() = default;
  YAML::Node findNode(const std::string& path, bool& found) const;
  void validateSchema() const;

  YAML::Node root_;
  std::string source_ = "<string>";
  std::string config_dir_;
};

}  // namespace frehg2

#endif  // FREHG2_IO_CONFIG_HPP
