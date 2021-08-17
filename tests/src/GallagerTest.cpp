#include <LDPCm/Gallager.hpp>
#include <catch2/catch.hpp>

using namespace LDPCm;

TEST_CASE("Gallager Gallager()", "[ldpc]") {
  SECTION("Bad (N, K) pairs") {
    // N > K
    CHECK_THROWS(Gallager(100, 100, 100, 5));
    // M multiple of 3
    CHECK_THROWS(Gallager(100, 20, 80, 3));
    // weightColumns * n / m must be an integer
    CHECK_THROWS(Gallager(100, 20, 80, 5));
  }
  SECTION("Matrix sizes initializations") {
    Gallager ldpc(100, 80, 20, 5);
    REQUIRE(ldpc.n == 100);
    REQUIRE(ldpc.k == 80);
    REQUIRE(ldpc.m == 20);
    REQUIRE(ldpc.weightColumns == 5);
    REQUIRE(ldpc.weightRows == 25);
    REQUIRE(ldpc.H.size() == 20);
    REQUIRE(ldpc.H[0].size() == 100);
    REQUIRE(ldpc.G.size() == 80);
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

TEST_CASE("Gallager generateH()", "[ldpc]") {
  SECTION("Small H") {
    Gallager ldpc(9, 6, 3, 3);
    ldpc.generateH();
    size_t ones = 0;
    for (auto row : ldpc.H)
      for (auto col : row) {
        ones += col;
      }
    REQUIRE(ones == ldpc.n * ldpc.weightColumns);
  }
  SECTION("Big H") {
    Gallager ldpc(100, 80, 20, 5);
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
  SECTION("Test numberOfRowsInH != m") {
    // Because Gallager was not good enough
    Gallager ldpc(100, 80, 20, 5);
    CHECK_THROWS(ldpc.generateHG());
  }
  Gallager ldpc(100, 80, 100, 5);
  SECTION("Test H") {
    size_t ones = 0;
    for (auto row : ldpc.H)
      for (auto col : row) {
        ones += col;
      }
    REQUIRE(ones == ldpc.n * ldpc.weightColumns);
  }
  SECTION("Test HRowEchelon") {
    for (size_t i = 0; i < ldpc.m; i++)
      for (size_t j = 0; j < ldpc.m; j++)
        if (i == j)
          REQUIRE(ldpc.HRowEchelon[i][j + ldpc.k] == 1);
        else
          REQUIRE(ldpc.HRowEchelon[i][j + ldpc.k] == 0);
  }
  SECTION("Test G") {
    for (size_t i = 0; i < ldpc.k; i++)
      for (size_t j = 0; j < ldpc.k; j++)
        if (i == j)
          REQUIRE(ldpc.G[i][j] == 1);
        else
          REQUIRE(ldpc.G[i][j] == 0);
  }
}
