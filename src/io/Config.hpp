#ifndef FREHG2_IO_CONFIG_HPP
#define FREHG2_IO_CONFIG_HPP

#include "core/types.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace frehg2 {

class Config {
public:
    explicit Config(const std::string& filename);

    bool has(const std::string& path) const;
    int indexedCount(const std::string& path) const;
    const std::string& filename() const noexcept;
    const std::string& baseDirectory() const noexcept;
    void validate() const;

    template <typename T>
    T get(const std::string& path) const;

    template <typename T>
    T getOr(const std::string& path, const T& default_value) const
    {
        return has(path) ? get<T>(path) : default_value;
    }

    GridSpec gridSpec() const;
    TimeSpec timeSpec() const;

private:
    std::string filename_;
    std::string base_directory_;
    std::unordered_map<std::string, std::string> values_;

    const std::string& requireRaw(const std::string& path) const;
};

template <>
std::string Config::get<std::string>(const std::string& path) const;

template <>
int Config::get<int>(const std::string& path) const;

template <>
real Config::get<real>(const std::string& path) const;

template <>
bool Config::get<bool>(const std::string& path) const;

template <>
std::vector<int> Config::get<std::vector<int>>(const std::string& path) const;

template <>
std::vector<real> Config::get<std::vector<real>>(const std::string& path) const;

template <>
std::vector<std::string> Config::get<std::vector<std::string>>(const std::string& path) const;

}  // namespace frehg2

#endif  // FREHG2_IO_CONFIG_HPP
