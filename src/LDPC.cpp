#include <LDPCm/LDPC.hpp>
#include <LDPCm/Uniform.hpp>
#include <stdexcept>
#include <algorithm>
#include <iostream>

using namespace LDPCm;
static Uniform uniform;

LDPC::LDPC(size_t _n, size_t _k, size_t _numberOfRowsInH)
    : n{_n},
      k{_k},
      m{_n - k},
      numberOfRowsInH{_numberOfRowsInH},
      H(numberOfRowsInH, std::vector<uint8_t>(n)),
      HRowEchelon(m, std::vector<uint8_t>(n)),
      G(k, std::vector<uint8_t>(n)),
      codewordsLinks(n, std::vector<size_t>()),
      paritiesLinks(numberOfRowsInH, std::vector<size_t>()),
      message(n),
      syndrome(numberOfRowsInH),
      codeword(n) {}

const std::vector<std::vector<uint8_t>> &LDPC::generateHRowEchelon() {
  auto HRowEchelonTemp = H;
  auto HTemp           = H;
  size_t i, j, r = 0, last = numberOfRowsInH;

  // iterate over rows
  while (true) {
    while (r < last) {
      // find the first "1" in this row after position r
      for (i = r; i < n; i++)
        if (HRowEchelonTemp[r][i])
          break;

      // row has no 1s after position r: exchange with another one
      if (i == n) {
        last--;

        // swap rows
        for (j = 0; j < n; j++)
          std::swap(HRowEchelonTemp[last][j], HRowEchelonTemp[r][j]);
      } else
        break;
    }

    // break if we reached the last row
    if (r == last)
      break;

    // swap column r and i to put a pivot in row r.
    for (j = 0; j < numberOfRowsInH; j++) {
      std::swap(HRowEchelonTemp[j][i], HRowEchelonTemp[j][r]);
      std::swap(H[j][i], H[j][r]);
    }

    // for all rows (except pivot)
    for (i = 0; i < numberOfRowsInH; i++)
      if ((i != r) && HRowEchelonTemp[i][r])
        for (j = r; j < n; j++)
          HRowEchelonTemp[i][j] = HRowEchelonTemp[i][j] ^ HRowEchelonTemp[r][j];

    r++;
  }

  // I|P to P|I
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      HRowEchelon[i][j] = HRowEchelonTemp[i][(j + m) % n];
      H[i][j]           = HTemp[i][(j + m) % n];
    }
  generateLinks();
  return HRowEchelon;
}

void LDPC::generateLinks() {
  for (size_t i = 0; i < numberOfRowsInH; i++)
    for (size_t j = 0; j < n; j++)
      if (H[i][j]) {
        codewordsLinks[j].push_back(i);
        paritiesLinks[i].push_back(j);
      }
}

const std::vector<std::vector<uint8_t>> &LDPC::generateG() {
  size_t i, j;
  for (i = 0; i < k; i++) {
    // I
    G[i][i] = 1;
    // P
    for (j = 0; j < n - k; j++)
      G[i][j + k] = HRowEchelon[j][i];
  }
  return G;
}

const std::vector<uint8_t> &LDPC::checkSyndrome(const std::vector<uint8_t> &_codeword) {
  size_t i, j;
  uint8_t temp;
  for (i = 0; i < numberOfRowsInH; i++) {
    temp = 0;
    for (j = 0; j < n; j++)
      if (H[i][j] && _codeword[j])
        temp ^= 1;
    syndrome[i] = temp;
  }
  return syndrome;
}

const std::vector<uint8_t> &LDPC::getCodeword(const std::vector<uint8_t> &_message) {
  size_t i, j, temp;

  for (i = 0; i < n; i++) {
    temp = 0;
    for (j = 0; j < k; j++)
      temp ^= _message[j] * G[j][i];
    codeword[i] = temp;
  }
  return codeword;
}

std::vector<uint8_t> LDPC::decoderBitFlip(std::vector<uint8_t> _message, std::vector<uint8_t> _syndrome, size_t _maxNumberOfIterations) {
  // std::vector<uint8_t> errorMessage(n);
  // std::vector<uint8_t> errorSyndrome(m);
  // std::vector<uint8_t> newSyndrome(_syndrome);
  // bool success;
  // size_t i, j, maxPositionMessage, maxPositionSyndrome;
  // while (0 < _maxNumberOfIterations--) {
  //   newSyndrome = getSyndrome(_message);
  //   success     = true;
  //   for (i = 0; i < m; i++)
  //     if (newSyndrome[i] ^ _syndrome[i]) {
  //       success = false;
  //       for (j = 0; j < n; j++)
  //         errorMessage[j] += H[i][j];
  //       errorSyndrome[i]++;
  //     }
  //   if (success)
  //     return _message;

  //   maxPositionMessage = 0;
  //   for (i = 0; i < n; i++) {
  //     if (errorMessage[i] > errorMessage[maxPositionMessage])
  //       maxPositionMessage = i;
  //     else
  //       errorMessage[i] = 0;
  //   }

  //   maxPositionSyndrome = 0;
  //   for (i = 0; i < m; i++) {
  //     if (errorSyndrome[i] > errorSyndrome[maxPositionSyndrome])
  //       maxPositionSyndrome = i;
  //     else
  //       errorSyndrome[i] = 0;
  //   }
  //   if (errorSyndrome[maxPositionSyndrome] > errorMessage[maxPositionMessage])
  //     _syndrome[maxPositionSyndrome] = _syndrome[maxPositionSyndrome] ? 0 : 1;
  //   else
  //     _message[maxPositionMessage] = _message[maxPositionMessage] ? 0 : 1;

  //   errorSyndrome[maxPositionSyndrome] = 0;
  //   errorMessage[maxPositionMessage]   = 0;
  // }
  return _message;
}

std::vector<uint8_t> LDPC::decoderBealivePropagation(std::vector<uint8_t> _codeword, size_t _maxNumberOfIterations) {
  std::vector<uint8_t> newSyndrome(numberOfRowsInH);
  std::vector<float> r(n);  // loglikelihood(0.005)
  codeword = _codeword;
  bool success;
  size_t i, j, ii, jj;
  float temp, errorProbability = 0.01;
  for (i = 0; i < n; i++)
    r[i] = codeword[i] ? loglikelihood(errorProbability) : -loglikelihood(errorProbability);

  std::vector<std::vector<float>> M(numberOfRowsInH, r);
  std::vector<std::vector<float>> E(numberOfRowsInH, std::vector<float>(n));
  std::vector<float> L(n);
  std::cout << "M matrix: " << M.size() << " " << M[0].size() << std::endl;
  std::cout << "E matrix: " << E.size() << " " << E[0].size() << std::endl;
  std::cout << "numberOfRowsInH: " << numberOfRowsInH << " n: " << n << std::endl;

  while (0 < _maxNumberOfIterations--) {
    // extrinsic probabilities
    for (j = 0; j < numberOfRowsInH; j++) {
      for (i = 0; i < paritiesLinks[j].size(); i++) {
        temp = 1;
        for (ii = 0; ii < paritiesLinks[j].size(); ii++)
          if (paritiesLinks[j][ii] != paritiesLinks[j][i])
            temp *= tanh(M[j][paritiesLinks[j][ii]] / 2);

        E[j][paritiesLinks[j][i]] = log((1 + temp) / (1 - temp));
      }
    }

    // loglikelihood
    for (i = 0; i < n; i++) {  // problem
      temp = 0;
      for (j = 0; j < codewordsLinks[i].size(); j++)
        temp += E[codewordsLinks[i][j]][i];

      L[i]        = r[i] + temp;
      codeword[i] = (L[i] > 0) ? 0 : 1;
    }

    // exit condition
    newSyndrome = checkSyndrome(codeword);
    success     = true;
    for (i = 0; i < numberOfRowsInH; i++)
      if (newSyndrome[i]) {
        success = false;
        break;
      }

    if (success)
      return codeword;

    // intrinsic probabilities
    for (i = 0; i < n; i++) {
      for (j = 0; j < codewordsLinks[i].size(); j++) {  // problem
        temp = 0;
        for (jj = 0; jj < codewordsLinks[i].size(); jj++)
          if (codewordsLinks[i][j] != codewordsLinks[i][jj])
            temp += E[codewordsLinks[i][jj]][i];

        M[codewordsLinks[i][j]][i] = temp + r[i];
      }
    }
  }
  return codeword;
}
