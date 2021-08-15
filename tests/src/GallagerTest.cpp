#include <LDPC/Gallager.hpp>
#include <catch2/catch.hpp>

using namespace LDPC;

TEST_CASE("Gallager Gallager()", "[ldpc]") {
  SECTION("Bad (N, K) pairs") {
    // N > K
    CHECK_THROWS(Gallager(100, 100, 5));
    // M multiple of 3
    CHECK_THROWS(Gallager(100, 20, 3));
  }
  SECTION("Matrix sizes initializations") {
    Gallager ldpc(100, 20, 5);
    REQUIRE(ldpc.H.size() == 80);
    REQUIRE(ldpc.H[0].size() == 100);
    REQUIRE(ldpc.G.size() == 20);
    REQUIRE(ldpc.G[0].size() == 100);
    size_t zero = 0;
    for (auto row : ldpc.H)
      for (auto col : row) {
        zero += col;
      }
    REQUIRE(zero == 0);
    zero = 0;
    for (auto row : ldpc.G)
      for (auto col : row) {
        zero += col;
      }
    REQUIRE(zero == 0);
  }
}

// TEST_CASE("Gallager setNumber(int _number)", "[ldpc]") {
//   Gallager ldpc;
//   ldpc.setNumber(1);
//   REQUIRE(ldpc.getNumber() == 1);
// }