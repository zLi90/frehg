// Minimal, self-contained SHA-256 (FIPS 180-4) for config provenance hashing (P3.2).
//
// Used to compute `config_sha256` (a stable hash of the resolved configuration string)
// embedded in every output file and checkpoint so a run is reproducible/traceable. This
// is a small dependency-free implementation; it is verified against the standard NIST
// test vectors in tests/io/test_sha256.cpp.
#ifndef FREHG2_IO_SHA256_HPP
#define FREHG2_IO_SHA256_HPP

#include <cstdint>
#include <string>

namespace frehg2 {

// Returns the lowercase hex SHA-256 digest of the input bytes.
std::string sha256Hex(const std::string& data);

}  // namespace frehg2

#endif  // FREHG2_IO_SHA256_HPP
