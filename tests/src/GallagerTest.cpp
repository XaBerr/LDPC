#include <LDPC/Gallager.hpp>
#include <catch2/catch.hpp>

using namespace LDPC;

TEST_CASE("Gallager Gallager()", "[ldpc]") {
  SECTION("Bad (N, K) pairs") {
    // N > K
    CHECK_THROWS(Gallager(100, 100, 5));
    // M multiple of 3
    CHECK_THROWS(Gallager(100, 20, 3));
    // weightColumns * n / m must be an integer
    CHECK_THROWS(Gallager(100, 20, 5));
  }
  SECTION("Matrix sizes initializations") {
    Gallager ldpc(100, 80, 5);
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

TEST_CASE("Gallager generateHRowEchelon()", "[ldpc]") {
  SECTION("Small H") {
    Gallager ldpc(8, 4, 2);
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
  SECTION("Big H") {
    Gallager ldpc(100, 80, 2);
    size_t i;
    for (i = 0; i < 100; i++) {
      ldpc.generateH();
      ldpc.generateHRowEchelon();
      if (ldpc.HRowEchelon[ldpc.m - 1][ldpc.n - 1])
        break;
      for (auto &row : ldpc.H)
        for (auto &col : row)
          col = 0;
    }
    // Error: cannot build a full linearly independent H.
    REQUIRE(i < 100);
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
}

// TEST_CASE("Gallager generateG()", "[ldpc]") {
//   Gallager ldpc(8, 4, 2);
//   ldpc.H = {
//       {1, 1, 1, 1, 0, 0, 0, 0},
//       {0, 0, 0, 0, 1, 1, 1, 1},

//       {0, 0, 1, 0, 1, 0, 1, 1},
//       {1, 1, 0, 0, 0, 1, 1, 0}};
//   ldpc.generateHRowEchelon();
//   ldpc.generateG();
//   for (size_t i = 0; i < ldpc.k; i++)
//     for (size_t j = 0; j < ldpc.k; j++)
//       if (i == j)
//         REQUIRE(ldpc.G[i][j] == 1);
//       else
//         REQUIRE(ldpc.G[i][j] == 0);
// }

// TEST_CASE("Gallager generateHG()", "[ldpc]") {
//   Gallager ldpc(100, 80, 5);
//   ldpc.generateHG();
//   SECTION("Test H") {
//     size_t ones = 0;
//     for (auto row : ldpc.H)
//       for (auto col : row) {
//         ones += col;
//       }
//     REQUIRE(ones == ldpc.n * ldpc.weightColumns);
//   }
//   SECTION("Test HRowEchelon") {
//     for (size_t i = 0; i < ldpc.m; i++)
//       for (size_t j = 0; j < ldpc.m; j++)
//         if (i == j)
//           REQUIRE(ldpc.HRowEchelon[i][j + ldpc.k] == 1);
//         else
//           REQUIRE(ldpc.HRowEchelon[i][j + ldpc.k] == 0);
//   }
//   SECTION("Test G") {
//     for (size_t i = 0; i < ldpc.k; i++)
//       for (size_t j = 0; j < ldpc.k; j++)
//         if (i == j)
//           REQUIRE(ldpc.G[i][j] == 1);
//         else
//           REQUIRE(ldpc.G[i][j] == 0);
//   }
// }

// TEST_CASE("Gallager getCodeword()", "[ldpc]") {
//   Gallager ldpc(8, 4, 2);
//   ldpc.H = {
//       {1, 1, 1, 1, 0, 0, 0, 0},
//       {0, 0, 0, 0, 1, 1, 1, 1},

//       {0, 0, 1, 0, 1, 0, 1, 1},
//       {1, 1, 0, 0, 0, 1, 1, 0}};
//   ldpc.generateHRowEchelon();
//   ldpc.generateG();

//   std::vector<uint8_t> message = {1, 1, 1, 1};
//   auto codeword                = ldpc.getCodeword(message);
//   auto syndrome                = ldpc.checkSyndrome(codeword);
//   size_t zero                  = 0;
//   for (auto pos : syndrome)
//     zero += pos;
//   REQUIRE(zero == 0);
// }
