#include <LDPCm/LDPC.hpp>
#include <catch2/catch.hpp>

using namespace LDPCm;

TEST_CASE("LDPC LDPC()", "[ldpc]") {
  SECTION("Bad (N, K) pairs") {
    // N > K
    CHECK_THROWS(LDPC(100, 100, 5));
    // M multiple of 3
    CHECK_THROWS(LDPC(100, 20, 3));
    // weightColumns * n / m must be an integer
    CHECK_THROWS(LDPC(100, 20, 5));
  }
  SECTION("Matrix sizes initializations") {
    LDPC ldpc(100, 80, 5);
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

TEST_CASE("LDPC generateHRowEchelon()", "[ldpc]") {
  LDPC ldpc(8, 4, 2);
  REQUIRE(ldpc.H.size() == 4);
  REQUIRE(ldpc.H[0].size() == 8);
  // rank 4 matrix checked with wolfram
  ldpc.H = {
      {1, 1, 1, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 1, 1, 1},

      {0, 0, 1, 0, 1, 0, 1, 1},
      {1, 1, 0, 0, 0, 1, 1, 0}};
  ldpc.generateHRowEchelon();
  for (size_t i = 0; i < ldpc.m; i++)
    for (size_t j = 0; j < ldpc.m; j++)
      if (i == j)
        REQUIRE(ldpc.HRowEchelon[i][j + ldpc.k] == 1);
      else
        REQUIRE(ldpc.HRowEchelon[i][j + ldpc.k] == 0);
  size_t ones = 0;
  for (auto row : ldpc.H)
    for (auto col : row) {
      ones += col;
    }
  REQUIRE(ones == ldpc.n * ldpc.weightColumns);
}

TEST_CASE("LDPC generateG()", "[ldpc]") {
  LDPC ldpc(8, 4, 2);
  ldpc.H = {
      {1, 1, 1, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 1, 1, 1},

      {0, 0, 1, 0, 1, 0, 1, 1},
      {1, 1, 0, 0, 0, 1, 1, 0}};
  ldpc.generateHRowEchelon();
  ldpc.generateG();
  for (size_t i = 0; i < ldpc.k; i++)
    for (size_t j = 0; j < ldpc.k; j++)
      if (i == j)
        REQUIRE(ldpc.G[i][j] == 1);
      else
        REQUIRE(ldpc.G[i][j] == 0);
}

TEST_CASE("LDPC getCodeword()", "[ldpc]") {
  LDPC ldpc(8, 4, 2);
  ldpc.H = {
      {1, 1, 1, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 1, 1, 1},

      {0, 0, 1, 0, 1, 0, 1, 1},
      {1, 1, 0, 0, 0, 1, 1, 0}};
  ldpc.generateHRowEchelon();
  ldpc.generateG();

  std::vector<uint8_t> message = {1, 1, 1, 1};
  auto codeword                = ldpc.getCodeword(message);
  auto syndrome                = ldpc.checkSyndrome(codeword);
  size_t zero                  = 0;
  for (auto pos : syndrome)
    zero += pos;
  REQUIRE(zero == 0);
}