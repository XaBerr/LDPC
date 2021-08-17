#include <LDPCm/LDPC.hpp>
#include <catch2/catch.hpp>

using namespace LDPCm;

TEST_CASE("LDPC LDPC()", "[ldpc]") {
  SECTION("Simple init") {
    LDPC ldpc(100, 80, 30);
    REQUIRE(ldpc.n == 100);
    REQUIRE(ldpc.k == 80);
    REQUIRE(ldpc.m == 20);
    REQUIRE(ldpc.numberOfRowsInH == 30);
    REQUIRE(ldpc.H.size() == 30);
    REQUIRE(ldpc.H[0].size() == 100);
    REQUIRE(ldpc.HRowEchelon.size() == 20);
    REQUIRE(ldpc.HRowEchelon[0].size() == 100);
    REQUIRE(ldpc.G.size() == 80);
    REQUIRE(ldpc.G[0].size() == 100);
  }
}

TEST_CASE("LDPC generateHRowEchelon()", "[ldpc]") {
  LDPC ldpc(8, 4, 4);
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
  REQUIRE(ones == 16);
}

TEST_CASE("LDPC generateG()", "[ldpc]") {
  LDPC ldpc(8, 4, 4);
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

TEST_CASE("LDPC getCodeword(message)", "[ldpc]") {
  LDPC ldpc(8, 4, 4);
  ldpc.H = {
      {1, 1, 1, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 1, 1, 1},

      {0, 0, 1, 0, 1, 0, 1, 1},
      {1, 1, 0, 0, 0, 1, 1, 0}};
  ldpc.generateHRowEchelon();
  ldpc.generateG();

  std::vector<uint8_t> message = {1, 1, 1, 1};
  auto codeword                = ldpc.getCodeword(message);
  std::vector<uint8_t> solution(8, 1);
  REQUIRE(codeword == solution);
}

TEST_CASE("LDPC checkSyndrome(codeword)", "[ldpc]") {
  LDPC ldpc(8, 4, 4);
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

TEST_CASE("LDPC decoderBealivePropagation(codeword,maxNumberOfIterations)", "[ldpc]") {
  LDPC ldpc(30, 22, 8);
  ldpc.H = {{1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0},
            {0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1},
            {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0},
            {0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},
            {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0},
            {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0},
            {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0},
            {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0}};
  ldpc.generateHRowEchelon();
  ldpc.generateG();

  std::vector<uint8_t> message{0, 0, 1};
  std::vector<uint8_t> codeword = ldpc.getCodeword(message);
  std::vector<uint8_t> codewordSolution{0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1};
  REQUIRE(codeword == codewordSolution);
  std::vector<uint8_t> syndrome = ldpc.checkSyndrome(codeword);
  std::vector<uint8_t> solution(8, 0);
  REQUIRE(syndrome == solution);
  codeword[1] = codeword[1] ? 0 : 1;
  codeword[3] = codeword[3] ? 0 : 1;
  codeword    = ldpc.decoderBealivePropagation(codeword);
  syndrome    = ldpc.checkSyndrome(codeword);
  REQUIRE(syndrome == solution);
}