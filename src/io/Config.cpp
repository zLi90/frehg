#include "io/Config.hpp"

#ifdef USE_YAML_CPP
#include <yaml-cpp/yaml.h>
#endif

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace frehg2 {

namespace {

std::string trim(const std::string& text)
{
    const auto first = std::find_if_not(text.begin(), text.end(), [](unsigned char ch) {
        return std::isspace(ch) != 0;
    });
    const auto last = std::find_if_not(text.rbegin(), text.rend(), [](unsigned char ch) {
        return std::isspace(ch) != 0;
    }).base();
    return first >= last ? std::string{} : std::string(first, last);
}

[[maybe_unused]] std::string stripComment(const std::string& line)
{
    bool in_quote = false;
    char quote = '\0';
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char ch = line[i];
        if ((ch == '"' || ch == '\'') && (i == 0 || line[i - 1] != '\\')) {
            if (!in_quote) {
                in_quote = true;
                quote = ch;
            } else if (quote == ch) {
                in_quote = false;
            }
        }
        if (ch == '#' && !in_quote) {
            return line.substr(0, i);
        }
    }
    return line;
}

std::string unquote(std::string value)
{
    value = trim(value);
    if (value.size() >= 2 &&
        ((value.front() == '"' && value.back() == '"') ||
         (value.front() == '\'' && value.back() == '\''))) {
        return value.substr(1, value.size() - 2);
    }
    return value;
}

std::vector<std::string> splitList(const std::string& raw)
{
    std::string value = trim(raw);
    if (value == "null" || value == "~") {
        return {};
    }
    if (value.size() >= 2 && value.front() == '[' && value.back() == ']') {
        value = value.substr(1, value.size() - 2);
    }

    std::vector<std::string> result;
    std::string item;
    bool in_quote = false;
    char quote = '\0';
    for (std::size_t i = 0; i < value.size(); ++i) {
        const char ch = value[i];
        if ((ch == '"' || ch == '\'') && (i == 0 || value[i - 1] != '\\')) {
            if (!in_quote) {
                in_quote = true;
                quote = ch;
            } else if (quote == ch) {
                in_quote = false;
            }
        }
        if (ch == ',' && !in_quote) {
            result.push_back(unquote(item));
            item.clear();
        } else {
            item.push_back(ch);
        }
    }
    if (!trim(item).empty()) {
        result.push_back(unquote(item));
    }
    return result;
}

[[maybe_unused]] int indentation(const std::string& line)
{
    int count = 0;
    for (const char ch : line) {
        if (ch != ' ') {
            break;
        }
        ++count;
    }
    return count;
}

[[maybe_unused]] std::string joinPath(
    const std::vector<std::string>& stack, int depth, const std::string& key)
{
    std::string path;
    for (int i = 0; i < depth; ++i) {
        if (!path.empty()) {
            path += ".";
        }
        path += stack[static_cast<std::size_t>(i)];
    }
    if (!path.empty()) {
        path += ".";
    }
    path += key;
    return path;
}

[[maybe_unused]] std::unordered_map<std::string, std::string> parseSimpleYaml(
    const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open configuration file: " + filename);
    }

    std::unordered_map<std::string, std::string> values;
    std::vector<std::string> stack;
    std::string line;
    while (std::getline(input, line)) {
        line = stripComment(line);
        if (trim(line).empty()) {
            continue;
        }

        const int depth = indentation(line) / 2;
        const std::string content = trim(line);
        if (content.rfind("- ", 0) == 0) {
            continue;
        }

        const auto colon = content.find(':');
        if (colon == std::string::npos) {
            throw std::runtime_error("invalid YAML line: " + content);
        }

        const std::string key = trim(content.substr(0, colon));
        const std::string value = trim(content.substr(colon + 1));

        if (static_cast<int>(stack.size()) <= depth) {
            stack.resize(static_cast<std::size_t>(depth + 1));
        }
        stack[static_cast<std::size_t>(depth)] = key;

        if (value.empty()) {
            values[joinPath(stack, depth, key)] = "__section__";
            continue;
        }

        values[joinPath(stack, depth, key)] = unquote(value);
    }

    return values;
}

#ifdef USE_YAML_CPP
std::string yamlScalarToString(const YAML::Node& node)
{
    if (!node || node.IsNull()) {
        return "null";
    }
    if (!node.IsScalar()) {
        return node.IsMap() ? "{}" : "[]";
    }
    return node.as<std::string>();
}

std::string yamlSequenceToString(const YAML::Node& node)
{
    std::string value = "[";
    for (std::size_t i = 0; i < node.size(); ++i) {
        if (i > 0) {
            value += ", ";
        }
        value += yamlScalarToString(node[i]);
    }
    value += "]";
    return value;
}

void flattenYamlNode(
    const std::string& path,
    const YAML::Node& node,
    std::unordered_map<std::string, std::string>& values)
{
    if (node.IsMap()) {
        if (!path.empty()) {
            values[path] = "__section__";
        }
        for (const auto& item : node) {
            const std::string key = item.first.as<std::string>();
            const std::string child_path = path.empty() ? key : path + "." + key;
            flattenYamlNode(child_path, item.second, values);
        }
        return;
    }

    if (node.IsSequence()) {
        values[path] = yamlSequenceToString(node);
        for (std::size_t i = 0; i < node.size(); ++i) {
            if (node[i].IsMap() || node[i].IsSequence()) {
                flattenYamlNode(path + "." + std::to_string(i), node[i], values);
            }
        }
        return;
    }

    values[path] = yamlScalarToString(node);
}

std::unordered_map<std::string, std::string> parseYamlCpp(const std::string& filename)
{
    std::unordered_map<std::string, std::string> values;
    flattenYamlNode("", YAML::LoadFile(filename), values);
    return values;
}
#endif

bool isNullLiteral(const std::string& value)
{
    return value == "null" || value == "~";
}

std::string parentDirectory(const std::string& filename)
{
    const auto separator = filename.find_last_of("/\\");
    if (separator == std::string::npos) {
        return ".";
    }
    if (separator == 0) {
        return filename.substr(0, 1);
    }
    return filename.substr(0, separator);
}

}  // namespace

Config::Config(const std::string& filename)
    : filename_(filename),
      base_directory_(parentDirectory(filename)),
      values_(
#ifdef USE_YAML_CPP
          parseYamlCpp(filename)
#else
          parseSimpleYaml(filename)
#endif
      )
{
}

bool Config::has(const std::string& path) const
{
    const auto iter = values_.find(path);
    return iter != values_.end() && !isNullLiteral(iter->second);
}

int Config::indexedCount(const std::string& path) const
{
    int count = 0;
    while (has(path + "." + std::to_string(count))) {
        ++count;
    }
    return count;
}

const std::string& Config::filename() const noexcept
{
    return filename_;
}

const std::string& Config::baseDirectory() const noexcept
{
    return base_directory_;
}

void Config::validate() const
{
    const std::vector<std::string> required_sections = {
        "simulation",
        "domain",
        "time",
        "surface_water",
        "groundwater",
        "solute",
        "output",
        "monitor",
    };
    for (const auto& section : required_sections) {
        if (values_.find(section) == values_.end()) {
            throw std::runtime_error("missing required YAML section: " + section);
        }
    }

    const std::vector<std::string> required_fields = {
        "simulation.id",
        "domain.nx",
        "domain.ny",
        "domain.dx",
        "domain.dy",
        "domain.dz",
        "domain.dz_incre",
        "domain.mpi_nx",
        "domain.mpi_ny",
        "time.dt",
        "time.Tend",
        "time.dt_out",
        "surface_water.enable",
        "surface_water.bc_type",
        "groundwater.enable",
        "groundwater.bc_type",
        "solute.enable",
        "output.filename",
        "monitor.n_monitor",
    };
    for (const auto& field : required_fields) {
        if (values_.find(field) == values_.end()) {
            throw std::runtime_error("missing required YAML field: " + field);
        }
    }

    if (get<std::vector<int>>("surface_water.bc_type").size() != 4) {
        throw std::runtime_error("surface_water.bc_type must have 4 entries");
    }
    if (get<std::vector<int>>("groundwater.bc_type").size() != 6) {
        throw std::runtime_error("groundwater.bc_type must have 6 entries");
    }
}

GridSpec Config::gridSpec() const
{
    GridSpec spec;
    spec.nx = get<int>("domain.nx");
    spec.ny = get<int>("domain.ny");
    spec.nz = getOr<int>("domain.nz", 1);
    if (spec.nz <= 0) {
        spec.nz = 1;
    }
    spec.dx = get<real>("domain.dx");
    spec.dy = get<real>("domain.dy");
    spec.dz = get<real>("domain.dz");
    spec.dz_multiplier = get<real>("domain.dz_incre");
    return spec;
}

TimeSpec Config::timeSpec() const
{
    TimeSpec spec;
    spec.start = 0.0;
    spec.end = get<real>("time.Tend");
    spec.dt = get<real>("time.dt");
    spec.output_interval = get<real>("time.dt_out");
    return spec;
}

const std::string& Config::requireRaw(const std::string& path) const
{
    const auto iter = values_.find(path);
    if (iter == values_.end() || isNullLiteral(iter->second)) {
        throw std::runtime_error("missing configuration value: " + path);
    }
    return iter->second;
}

template <>
std::string Config::get<std::string>(const std::string& path) const
{
    return unquote(requireRaw(path));
}

template <>
int Config::get<int>(const std::string& path) const
{
    return std::stoi(requireRaw(path));
}

template <>
real Config::get<real>(const std::string& path) const
{
    return std::stod(requireRaw(path));
}

template <>
bool Config::get<bool>(const std::string& path) const
{
    const std::string value = requireRaw(path);
    if (value == "true" || value == "True" || value == "1") {
        return true;
    }
    if (value == "false" || value == "False" || value == "0") {
        return false;
    }
    throw std::runtime_error("invalid boolean configuration value: " + path);
}

template <>
std::vector<int> Config::get<std::vector<int>>(const std::string& path) const
{
    std::vector<int> result;
    for (const auto& item : splitList(requireRaw(path))) {
        result.push_back(std::stoi(item));
    }
    return result;
}

template <>
std::vector<real> Config::get<std::vector<real>>(const std::string& path) const
{
    std::vector<real> result;
    for (const auto& item : splitList(requireRaw(path))) {
        result.push_back(std::stod(item));
    }
    return result;
}

template <>
std::vector<std::string> Config::get<std::vector<std::string>>(const std::string& path) const
{
    return splitList(requireRaw(path));
}

}  // namespace frehg2
