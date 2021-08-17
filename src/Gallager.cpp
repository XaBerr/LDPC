#include <LDPCm/Gallager.hpp>
#include <LDPCm/Uniform.hpp>
#include <stdexcept>
#include <algorithm>
#include <iostream>

using namespace LDPCm;
static Uniform uniform;

Gallager::Gallager(size_t _n, size_t _k, size_t _weightColumns)
    : LDPC(_n, _k, _weightColumns) {}

const std::vector<std::vector<uint8_t>> &Gallager::generateH() {
  size_t i, j, u, k, v;
  std::vector<std::vector<size_t>> shuffles(weightColumns, std::vector<size_t>(n));
  for (u = 0; u < weightColumns; u++) {
    for (i = 0; i < n; ++i) {
      shuffles[u][i] = i;
    }
    if (u > 0)
      std::shuffle(shuffles[u].begin(), shuffles[u].end(), uniform.getGenerator());
  }
  for (i = 0; i < m / weightColumns; i++)
    for (j = 0; j < weightColumns * n / m; j++) {
      k = j + i * weightColumns * n / m;
      for (u = 0; u < weightColumns; u++) {
        v = i + u * m / weightColumns;
        // bands
        H[v][shuffles[u][k]] = 1;
        codewordsLinks[shuffles[u][k]].push_back(v);
        paritiesLinks[v].push_back(shuffles[u][k]);
      }
    }
  return H;
}

void Gallager::generateHG(int _maxNumberOfIterations) {
  size_t i;
  for (i = 0; i < _maxNumberOfIterations; i++) {
    generateH();
    generateHRowEchelon();
    if (HRowEchelon[m - 1][n - 1])
      break;
    for (auto &row : H)
      for (auto &col : row)
        col = 0;
  }
  if (i >= 100)
    throw std::runtime_error("Error: cannot build a full linearly independent H.");
  generateG();
}
