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

std::vector<int> Gallager::decoderBitFlip(std::vector<int> _message, std::vector<int> _syndrome, int _maxNumberOfIterations) {
  std::vector<int> errorMessage(columns);
  std::vector<int> errorSyndrome(rows);
  std::vector<int> newSyndrome(_syndrome);
  bool success;
  int i, j, maxPositionMessage, maxPositionSyndrome;
  while (0 < _maxNumberOfIterations--) {
    newSyndrome = getSyndrome(_message);
    success     = true;
    for (i = 0; i < rows; i++)
      if (newSyndrome[i] ^ _syndrome[i]) {
        success = false;
        for (j = 0; j < columns; j++)
          errorMessage[j] += H[i][j];
        errorSyndrome[i]++;
      }
    if (success)
      return _message;

    maxPositionMessage = 0;
    for (i = 0; i < columns; i++) {
      if (errorMessage[i] > errorMessage[maxPositionMessage])
        maxPositionMessage = i;
      else
        errorMessage[i] = 0;
    }

    maxPositionSyndrome = 0;
    for (i = 0; i < rows; i++) {
      if (errorSyndrome[i] > errorSyndrome[maxPositionSyndrome])
        maxPositionSyndrome = i;
      else
        errorSyndrome[i] = 0;
    }
    if (errorSyndrome[maxPositionSyndrome] > errorMessage[maxPositionMessage])
      _syndrome[maxPositionSyndrome] = _syndrome[maxPositionSyndrome] ? 0 : 1;
    else
      _message[maxPositionMessage] = _message[maxPositionMessage] ? 0 : 1;

    errorSyndrome[maxPositionSyndrome] = 0;
    errorMessage[maxPositionMessage]   = 0;
  }
  return _message;
}

std::vector<int> Gallager::decoderBealivePropagation(std::vector<int> _message, std::vector<int> _syndrome, int _maxNumberOfIterations) {
  std::vector<float> r(columns, loglikelihood(0.5));
  std::vector<std::vector<float>> M(rows, std::vector<float>(columns, loglikelihood(0.5)));
  std::vector<std::vector<float>> E(rows, std::vector<float>(columns, 0));
  std::vector<float> L(columns);
  bool success;
  int i, j, jj;
  float temp;
  while (0 < _maxNumberOfIterations--) {
    // extrinsic probabilities
    for (j = 0; j < rows; j++) {
      temp = 1;
      for (i = 0; i < columns; i++) {
        temp *= tanh(M[j][i] / 2);
      }
      E[j][i] = log((1 + temp) / (1 - temp));
    }
    // loglikelihood
    success = true;
    for (i = 0; i < columns; i++) {
      temp = 0;
      for (j = 0; j < rows; j++)
        temp += E[j][i];

      L[i] = r[i] + temp;
      if (L[i] > 0 != _message[i])
        success = false;
      _message[i] = (L[i] > 0) ? 0 : 1;
    }
    if (success)
      return _message;
    // intrinsic probabilities
    for (i = 0; i < columns; i++) {
      for (j = 0; j < rows; j++) {
        temp = 0;
        for (jj = 0; jj < rows; jj++)
          if (j != jj)
            temp += E[jj][i];

        M[j][i] = temp + r[i];
      }
    }
  }
  return _message;
}