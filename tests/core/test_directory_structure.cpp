#include <filesystem>
#include <string>
#include <vector>

int main()
{
    const std::filesystem::path root = FREHG2_SOURCE_DIR;
    const std::vector<std::string> required_directories = {
        "benchmarks/b1-sw",
        "benchmarks/b2-gw",
        "cmake",
        "docs/legacy_audit",
        "external",
        "scripts",
        "src/core",
        "src/io",
        "src/swe",
        "src/re",
        "src/solute",
        "src/coupling",
        "src/bc",
        "tests/core",
        "tests/io",
        "tests/swe",
        "tests/re",
        "tests/solute",
        "tests/coupling",
        "tests/integration",
    };

    for (const auto& directory : required_directories) {
        if (!std::filesystem::is_directory(root / directory)) {
            return 1;
        }
    }

    return 0;
}
