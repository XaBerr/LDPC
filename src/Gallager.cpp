#include <LDPC/Gallager.hpp>
#include <LDPC/Uniform.hpp>
#include <stdexcept>
#include <algorithm>
#include <iostream>

using namespace LDPC;
static Uniform uniform;

Gallager::Gallager(int _n, int _k)
    : n{validateN(_n, _k)},
      k{validateK(_n, _k)},
      rows{_n - _k},
      columns{_k},
      H(rows, std::vector<int>(columns)),
      message(columns),
      syndrome(rows),
      codeword(_n) {
  std::vector<int> shuffle1(columns);
  std::vector<int> shuffle2(columns);
  for (int i = 0; i < columns; ++i) {
    shuffle1[i] = i;
    shuffle2[i] = i;
  }
  std::shuffle(shuffle1.begin(), shuffle1.end(), uniform.getGenerator());
  std::shuffle(shuffle2.begin(), shuffle2.end(), uniform.getGenerator());
  for (size_t i = 0; i < rows / 3; i++)
    for (size_t j = 0; j < 3 * columns / rows; j++) {
      // first band
      H[i][j + i * 3 * columns / rows] = 1;
      // second band
      H[i + 1 * rows / 3][shuffle1[j + i * 3 * columns / rows]] = 1;
      // third band
      H[i + 2 * rows / 3][shuffle2[j + i * 3 * columns / rows]] = 1;
    }
}

int Gallager::validateN(int _n, int _k) {
  if (_n <= 0)
    throw std::invalid_argument("Error: n must be > 0.");
  return _n;
}

int Gallager::validateK(int _n, int _k) {
  if (_k <= 0)
    throw std::invalid_argument("Error: k must be > 0.");
  if (_n <= _k)
    throw std::invalid_argument("Error: n must be > k.");
  int m = _n - _k;
  if (m % 3)
    throw std::invalid_argument("Error: n - k must be multiple of 3.");
  if ((3 * _n / m) % 1)
    throw std::invalid_argument("Error: 3 * n / (n -k) must an integer.");
  return _k;
}

const std::vector<int> &Gallager::getSyndrome(const std::vector<int> &_message) {
  int temp;
  for (size_t i = 0; i < rows; i++) {
    temp = 0;
    for (size_t j = 0; j < columns; j++)
      temp ^= H[i][j] * _message[j];
    syndrome[i] = temp;
  }
  return syndrome;
}

const std::vector<int> &Gallager::decode(const std::vector<int> &_message, const std::vector<int> &_syndrome) {
  return _message;
}