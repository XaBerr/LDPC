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
      m{_n - _k},
      H(m, std::vector<uint8_t>(n)),
      G(m, std::vector<uint8_t>(n)),
      codewordsLinks(n, std::vector<uint8_t>()),
      paritiesLinks(m, std::vector<uint8_t>()),
      message(n),
      syndrome(m),
      codeword(n) {
  std::vector<uint8_t> shuffle1(n);
  std::vector<uint8_t> shuffle2(n);
  for (int i = 0; i < n; ++i) {
    shuffle1[i] = i;
    shuffle2[i] = i;
  }
  std::shuffle(shuffle1.begin(), shuffle1.end(), uniform.getGenerator());
  std::shuffle(shuffle2.begin(), shuffle2.end(), uniform.getGenerator());
  int k;
  for (size_t i = 0; i < m / 3; i++)
    for (size_t j = 0; j < 3 * n / m; j++) {
      k = j + i * 3 * n / m;
      // first band
      H[i][k] = 1;
      codewordsLinks[k].push_back(i);
      paritiesLinks[i].push_back(k);
      // second band
      H[i + 1 * m / 3][shuffle1[k]] = 1;
      codewordsLinks[shuffle1[k]].push_back(i + 1 * m / 3);
      paritiesLinks[i + 1 * m / 3].push_back(shuffle1[k]);
      // third band
      H[i + 2 * m / 3][shuffle2[k]] = 1;
      codewordsLinks[shuffle2[k]].push_back(i + 2 * m / 3);
      paritiesLinks[i + 2 * m / 3].push_back(shuffle2[k]);
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
  if ((_n - _k) % 3)
    throw std::invalid_argument("Error: n - k must be multiple of 3.");
  if ((3 * _n / (_n - _k)) % 1)
    throw std::invalid_argument("Error: 3 * n / (n -k) must an integer.");
  return _k;
}

const std::vector<uint8_t> &Gallager::getSyndrome(const std::vector<uint8_t> &_message) {
  int temp;
  for (int i = 0; i < m; i++) {
    temp = 0;
    // I use only the first part of H
    for (int j = 0; j < k; j++)
      temp ^= H[i][j] * _message[j];
    syndrome[i] = temp;
  }
  return syndrome;
}

std::vector<uint8_t> Gallager::decoderBitFlip(std::vector<uint8_t> _message, std::vector<uint8_t> _syndrome, int _maxNumberOfIterations) {
  // std::vector<uint8_t> errorMessage(n);
  // std::vector<uint8_t> errorSyndrome(m);
  // std::vector<uint8_t> newSyndrome(_syndrome);
  // bool success;
  // int i, j, maxPositionMessage, maxPositionSyndrome;
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

std::vector<uint8_t> Gallager::decoderBealivePropagation(std::vector<uint8_t> _message, std::vector<uint8_t> _syndrome, int _maxNumberOfIterations) {
  std::vector<uint8_t> newSyndrome(_syndrome);
  std::vector<uint8_t> codeword(n);
  // std::vector<float> r(n);  // loglikelihood(0.005)
  // bool success;
  // int i, j, ii, jj;
  // float temp, errorProbability = 0.01;
  // for (i = 0; i < n; i++) {
  //   codeword[i] = _message[i];
  //   r[i]        = _message[i] ? loglikelihood(errorProbability) : -loglikelihood(errorProbability);
  // }
  // for (j = 0; j < m; j++) {
  //   codeword[j + n] = _message[j];
  //   r[j + n]        = _syndrome[j] ? loglikelihood(errorProbability) : -loglikelihood(errorProbability);
  // }

  // std::vector<std::vector<float>> M(m, r);
  // std::vector<std::vector<float>> E(m, std::vector<float>(n));
  // std::vector<float> L(n);
  // std::cout << "M matrix: " << M.size() << " " << M[0].size() << std::endl;
  // std::cout << "E matrix: " << E.size() << " " << E[0].size() << std::endl;
  // std::cout << "m: " << m << " n: " << n << std::endl;

  // while (0 < _maxNumberOfIterations--) {
  //   // extrinsic probabilities
  //   for (j = 0; j < m; j++) {
  //     for (i = 0; i < paritiesLinks[j].size(); i++) {
  //       temp = 1;
  //       for (ii = 0; ii < paritiesLinks[j].size(); ii++)
  //         if (paritiesLinks[j][ii] != paritiesLinks[j][i])
  //           temp *= tanh(M[j][paritiesLinks[j][ii]] / 2);

  //       E[j][paritiesLinks[j][i]] = log((1 + temp) / (1 - temp));
  //     }
  //   }

  //   // loglikelihood
  //   for (i = 0; i < n; i++) {  // problem
  //     temp = 0;
  //     for (j = 0; j < codewordsLinks[i].size(); j++)
  //       temp += E[codewordsLinks[i][j]][i];

  //     L[i]        = r[i] + temp;
  //     codeword[i] = (L[i] > 0) ? 0 : 1;
  //   }

  //   // exit condition
  //   newSyndrome = getSyndrome(codeword);
  //   success     = true;
  //   for (i = 0; i < m; i++)
  //     if (newSyndrome[i] == codeword[i + n]) {
  //       success = false;
  //       break;
  //     }

  //   if (success)
  //     return _message;

  //   // intrinsic probabilities
  //   for (i = 0; i < n; i++) {
  //     for (j = 0; j < codewordsLinks[i].size(); j++) {  // problem
  //       temp = 0;
  //       for (jj = 0; jj < codewordsLinks[i].size(); jj++)
  //         if (codewordsLinks[i][j] != codewordsLinks[i][jj])
  //           temp += E[codewordsLinks[i][jj]][i];

  //       M[codewordsLinks[i][j]][i] = temp + r[i];
  //     }
  //   }
  // }
  return codeword;
}

const std::vector<std::vector<uint8_t>> &Gallager::gaussianElimination() {
  H2 = H;
  int i, j, r = 0, last = n - k;

  // iterate over rows
  while (true) {
    while (r < last) {
      // find the first "1" in this row after position r
      for (i = r; i < n; i++)
        if (H2[r][i])
          break;

      // row has no 1s after position r: exchange with another one
      if (i == n) {
        last--;

        // swap rows
        for (j = 0; j < n; j++)
          std::swap(H2[last][j], H2[r][j]);
      } else
        break;
    }

    // break if we reached the last row
    if (r == last)
      break;

    // swap column r and i to put a pivot in row r.
    for (j = 0; j < n - k; j++) {
      std::swap(H2[j][i], H2[j][r]);
      std::swap(H[j][i], H[j][r]);
    }

    // for all rows (except pivot)
    for (i = 0; i < n - k; i++)
      if ((i != r) && H2[i][r])
        for (j = r; j < n; j++)
          H2[i][j] = H2[i][j] ^ H2[r][j];

    r++;
  }
  return H2;
}

const std::vector<std::vector<uint8_t>> &Gallager::generateG() {
  for (int i = 0; i < k; i++) {
    // P
    for (int j = 0; j < n - k; j++)
      G[i][j] = H2[j][i + (n - k)];

    // I
    for (int j = n - k; j < n; j++)
      if (j - (n - k) == i)
        G[i][j] = 1;
  }
  return G;
}
