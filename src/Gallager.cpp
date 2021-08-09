#include <LDPC/Gallager.hpp>
#include <LDPC/Uniform.hpp>
#include <stdexcept>

using namespace LDPC;
static Uniform uniform;

Gallager::Gallager(int _n, int _k)
    : n{validateN(_n, _k)},
      k{validateK(_n, _k)},
      rows{_k},
      columns{_n - _k},
      H(rows, std::vector<int>(columns, 0)),
      message{_n - _k, 0},
      syndrome{_k, 0},
      codeword{_n, 0} {
  // first band
  int m = _n - _k;
  for (size_t i = 0; i < rows / 3; i++)
    for (size_t j = 0; j < 3 * columns / m; j++)
      H[i][j + i * 3 * columns / m] = 1;
  // second band
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
    for (size_t j = 0; j < columns; j++) {
      temp += H[i][j] * _message[j];
    }
    syndrome[i] = temp;
  }

  return syndrome;
}

const std::vector<int> &Gallager::decode(const std::vector<int> &_message, const std::vector<int> &_syndrome) {
  return _message;
}