#include <LDPC/Gallager.hpp>
#include <stdexcept>

using namespace LDPC;

Gallager::Gallager(int _n, int _k)
    : n{validateN(_n, _k)},
      k{validateK(_n, _k)},
      rows{_k},
      columns{_n - _k},
      H(rows, std::vector<int>(columns, 0)),
      message{_n - _k, 0},
      syndrome{_k, 0},
      codeword{_n, 0} {
  if (_n <= _k)
    throw std::invalid_argument("Error: n must be > k.");
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