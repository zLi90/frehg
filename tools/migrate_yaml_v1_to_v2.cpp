// migrate_yaml_v1_to_v2 (P17, Task 17.3.4)
//
// Reads a v1/experimental YAML config, applies the frozen-v2 rename map (see
// src/io/YamlMigration.hpp / docs/yaml_schema_v2.md), and writes a v2 YAML.
//
// Usage:
//   migrate_yaml_v1_to_v2 <input.yaml> [output.yaml]
// With no output path, the v2 YAML is written to stdout.
#include <yaml-cpp/yaml.h>

#include <fstream>
#include <iostream>
#include <string>

#include "io/YamlMigration.hpp"

int main(int argc, char** argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "usage: " << argv[0] << " <input.yaml> [output.yaml]\n";
    return 2;
  }
  const std::string in_path = argv[1];

  YAML::Node in;
  try {
    in = YAML::LoadFile(in_path);
  } catch (const YAML::Exception& e) {
    std::cerr << "migrate_yaml_v1_to_v2: failed to parse '" << in_path << "': " << e.what() << "\n";
    return 1;
  }

  YAML::Node out;
  try {
    out = frehg2::migrateV1ToV2(in);
  } catch (const std::exception& e) {
    std::cerr << "migrate_yaml_v1_to_v2: migration failed: " << e.what() << "\n";
    return 1;
  }

  std::stringstream ss;
  ss << out << "\n";

  if (argc == 3) {
    const std::string out_path = argv[2];
    std::ofstream os(out_path);
    if (!os) {
      std::cerr << "migrate_yaml_v1_to_v2: cannot open output '" << out_path << "'\n";
      return 1;
    }
    os << ss.str();
    std::cerr << "migrate_yaml_v1_to_v2: wrote " << out_path << " (schema_version 2.0)\n";
  } else {
    std::cout << ss.str();
  }
  return 0;
}
