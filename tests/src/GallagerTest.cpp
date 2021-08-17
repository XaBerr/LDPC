#include <LDPCm/Gallager.hpp>
#include <catch2/catch.hpp>

using namespace LDPCm;

TEST_CASE("Gallager Gallager()", "[ldpc]") {
  Gallager ldpc(100, 80, 5);
  REQUIRE(true);
}

TEST_CASE("Gallager generateH()", "[ldpc]") {
  SECTION("Small H") {
    Gallager ldpc(9, 6, 3);
    ldpc.generateH();
    size_t ones = 0;
    for (auto row : ldpc.H)
      for (auto col : row) {
        ones += col;
      }
    REQUIRE(ones == ldpc.n * ldpc.weightColumns);
  }
  SECTION("Big H") {
    Gallager ldpc(100, 80, 5);
    ldpc.generateH();
    SECTION("Test H") {
      size_t ones = 0;
      for (auto row : ldpc.H)
        for (auto col : row) {
          ones += col;
        }
      REQUIRE(ones == ldpc.n * ldpc.weightColumns);
    }
  }
}

TEST_CASE("Gallager generateHG()", "[ldpc]") {
  Gallager ldpc(100, 80, 5);
  // Because Gallager was not good enough
  CHECK_THROWS(ldpc.generateHG());

  // SECTION("Test H") {
  //   size_t ones = 0;
  //   for (auto row : ldpc.H)
  //     for (auto col : row) {
  //       ones += col;
  //     }
  //   REQUIRE(ones == ldpc.n * ldpc.weightColumns);
  // }
  // SECTION("Test HRowEchelon") {
  //   for (size_t i = 0; i < ldpc.m; i++)
  //     for (size_t j = 0; j < ldpc.m; j++)
  //       if (i == j)
  //         REQUIRE(ldpc.HRowEchelon[i][j + ldpc.k] == 1);
  //       else
  //         REQUIRE(ldpc.HRowEchelon[i][j + ldpc.k] == 0);
  // }
  // SECTION("Test G") {
  //   for (size_t i = 0; i < ldpc.k; i++)
  //     for (size_t j = 0; j < ldpc.k; j++)
  //       if (i == j)
  //         REQUIRE(ldpc.G[i][j] == 1);
  //       else
  //         REQUIRE(ldpc.G[i][j] == 0);
  // }
}
