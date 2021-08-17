#include <LDPCm/Gallager.hpp>
#include <LDPCm/Uniform.hpp>
#include <stdexcept>
#include <algorithm>
#include <iostream>

using namespace LDPCm;
static Uniform uniform;

Gallager::Gallager(size_t _n, size_t _k, size_t _numberOfRowsInH, size_t _weightColumns)
    : LDPC(validateN(_n), validateK(_n, _k), _numberOfRowsInH),
      weightColumns{validateWeightColumns(_n, _k, _numberOfRowsInH, _weightColumns)},
      weightRows{weightColumns * n / numberOfRowsInH} {}

size_t Gallager::validateN(size_t _n) {
  if (_n <= 0)
    throw std::invalid_argument("Error: n must be > 0.");
  return _n;
}

size_t Gallager::validateK(size_t _n, size_t _k) {
  if (_k <= 0)
    throw std::invalid_argument("Error: k must be > 0.");
  if (_n <= _k)
    throw std::invalid_argument("Error: k must be < n.");
  return _k;
}

size_t Gallager::validateWeightColumns(size_t _n, size_t _k, size_t _numberOfRowsInH, size_t _weightColumns) {
  if ((_n - _k) > _numberOfRowsInH)
    throw std::invalid_argument("Error: (n - k) must be > that numberOfRowsInH.");
  if ((_numberOfRowsInH) % _weightColumns)
    throw std::invalid_argument("Error: numberOfRowsInH must be multiple of weightColumns.");
  float x = ((float)_weightColumns * _n / (_numberOfRowsInH));
  if (std::floor(x) != x)
    throw std::invalid_argument("Error: weightColumns * n / numberOfRowsInH must an integer.");
  return _weightColumns;
}

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
  for (i = 0; i < numberOfRowsInH / weightColumns; i++)
    for (j = 0; j < weightColumns * n / numberOfRowsInH; j++) {
      k = j + i * weightColumns * n / numberOfRowsInH;
      for (u = 0; u < weightColumns; u++) {
        v = i + u * numberOfRowsInH / weightColumns;
        // bands
        H[v][shuffles[u][k]] = 1;
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
