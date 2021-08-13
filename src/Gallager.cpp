#include <LDPC/Gallager.hpp>
#include <LDPC/Uniform.hpp>
#include <stdexcept>
#include <algorithm>
#include <iostream>

using namespace LDPC;
static Uniform uniform;

Gallager::Gallager(size_t _n, size_t _k, size_t _weightColumns)
    : n{validateN(_n)},
      k{validateK(_n, _k)},
      m{validateM(_n, _k, _weightColumns)},
      weightColumns{_weightColumns},
      weightRows{weightColumns * n / m},
      H(m, std::vector<uint8_t>(n)),
      HRowEchelon(),
      G(k, std::vector<uint8_t>(n)),
      codewordsLinks(n, std::vector<size_t>()),
      paritiesLinks(m, std::vector<size_t>()),
      message(n),
      syndrome(m),
      codeword(n) {
  size_t i;
  for (i = 0; i < 100; i++) {
    generateH();
    generateHRowEchelon();
    if (HRowEchelon[m - 1][m - 1])
      break;
  }
  if (i >= 100)
    throw std::runtime_error("Error: cannot build a full linearly independent H.");
  generateG();
}

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

size_t Gallager::validateM(size_t _n, size_t _k, size_t _weightColumns) {
  if ((_n - _k) % _weightColumns)
    throw std::invalid_argument("Error: n - k (i.e. m) must be multiple of weightColumns.");
  if ((_weightColumns * _n / (_n - _k)) % 1)
    throw std::invalid_argument("Error: weightColumns * n / (n -k) must an integer.");
  return _n - _k;
}

const std::vector<uint8_t> &Gallager::checkSyndrome(const std::vector<uint8_t> &_codeword) {
  size_t i, j, temp;
  for (i = 0; i < k; i++) {
    temp = 0;
    for (j = 0; j < n; j++)
      temp ^= H[i][j] * _codeword[j];
    syndrome[i] = temp;
  }
  return syndrome;
}

const std::vector<uint8_t> &Gallager::getCodeword(const std::vector<uint8_t> &_message) {
  size_t i, j, temp;

  for (i = 0; i < n; i++) {
    temp = 0;
    for (j = 0; j < k; j++)
      temp ^= _message[j] * G[j][i];
    codeword[i] = temp;
  }
  return codeword;
}

std::vector<uint8_t> Gallager::decoderBitFlip(std::vector<uint8_t> _message, std::vector<uint8_t> _syndrome, size_t _maxNumberOfIterations) {
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

std::vector<uint8_t> Gallager::decoderBealivePropagation(std::vector<uint8_t> _codeword, size_t _maxNumberOfIterations) {
  std::vector<uint8_t> newSyndrome(m);
  std::vector<float> r(n);  // loglikelihood(0.005)
  codeword = _codeword;
  bool success;
  size_t i, j, ii, jj;
  float temp, errorProbability = 0.01;
  for (i = 0; i < n; i++)
    r[i] = codeword[i] ? loglikelihood(errorProbability) : -loglikelihood(errorProbability);

  std::vector<std::vector<float>> M(m, r);
  std::vector<std::vector<float>> E(m, std::vector<float>(n));
  std::vector<float> L(n);
  std::cout << "M matrix: " << M.size() << " " << M[0].size() << std::endl;
  std::cout << "E matrix: " << E.size() << " " << E[0].size() << std::endl;
  std::cout << "m: " << m << " n: " << n << std::endl;

  while (0 < _maxNumberOfIterations--) {
    // extrinsic probabilities
    for (j = 0; j < m; j++) {
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
    for (i = 0; i < m; i++)
      if (!newSyndrome[i]) {
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

const std::vector<std::vector<uint8_t>> &Gallager::generateHRowEchelon() {
  HRowEchelon = H;
  size_t i, j, r = 0, last = n - k;

  // iterate over rows
  while (true) {
    while (r < last) {
      // find the first "1" in this row after position r
      for (i = r; i < n; i++)
        if (HRowEchelon[r][i])
          break;

      // row has no 1s after position r: exchange with another one
      if (i == n) {
        last--;

        // swap rows
        for (j = 0; j < n; j++)
          std::swap(HRowEchelon[last][j], HRowEchelon[r][j]);
      } else
        break;
    }

    // break if we reached the last row
    if (r == last)
      break;

    // swap column r and i to put a pivot in row r.
    for (j = 0; j < n - k; j++) {
      std::swap(HRowEchelon[j][i], HRowEchelon[j][r]);
      std::swap(H[j][i], H[j][r]);
    }

    // for all rows (except pivot)
    for (i = 0; i < n - k; i++)
      if ((i != r) && HRowEchelon[i][r])
        for (j = r; j < n; j++)
          HRowEchelon[i][j] = HRowEchelon[i][j] ^ HRowEchelon[r][j];

    r++;
  }
  return HRowEchelon;
}

const std::vector<std::vector<uint8_t>> &Gallager::generateG() {
  size_t i, j;
  for (i = 0; i < k; i++) {
    // I
    G[i][i] = 1;
    // P
    for (j = 0; j < n - k; j++)
      G[i][j + k] = HRowEchelon[j][i + (n - k)];
  }
  return G;
}
