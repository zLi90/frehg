// P3.2 provenance: SHA-256 used for config_sha256. Verify against NIST FIPS-180-4 vectors.
#define FREHG2_TEST_IMPL
#include "frehg2_test.hpp"
#include "io/Sha256.hpp"

using namespace frehg2;

TEST_CASE("sha256 known vectors") {
  REQUIRE(sha256Hex("") ==
          std::string("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"));
  REQUIRE(sha256Hex("abc") ==
          std::string("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"));
  REQUIRE(sha256Hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ==
          std::string("248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"));
}

TEST_CASE("sha256 is deterministic and sensitive to changes") {
  const std::string a = sha256Hex("frehg2 config v1");
  const std::string b = sha256Hex("frehg2 config v1");
  const std::string c = sha256Hex("frehg2 config v2");
  REQUIRE(a == b);
  REQUIRE(a != c);
  REQUIRE(a.size() == 64u);
}
