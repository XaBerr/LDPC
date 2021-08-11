#ifndef LDPC_GALLAGER_H
#define LDPC_GALLAGER_H

#include <vector>
#include <math.h>

namespace LDPC {
class Gallager {
  const int n;
  const int k;
  const int rows;
  const int columns;
  std::vector<std::vector<int>> H;
  std::vector<int> message;
  std::vector<int> syndrome;
  std::vector<int> codeword;

  int validateN(int _n, int _k);
  int validateK(int _n, int _k);
  static inline float loglikelihood(float d) { return log(d / (1 - d)); }
  // static inline float atanh(float f) { return log((1 + f) / (1 - f)); }
  static inline float phy(float x) { return log((exp(x) + 1) / (exp(x) - 1)); }

 public:
  Gallager(int _n, int _k);
  // const std::vector<int> &decode(const std::vector<int> &_codeword);
  // const std::vector<int> &encode(const std::vector<int> &_message);
  const std::vector<int> &getSyndrome(const std::vector<int> &_message);
  // const std::vector<int> &decode(const std::vector<int> &_message, const std::vector<int> &_syndrome);
  std::vector<int> decoderBitFlip(std::vector<int> _message, std::vector<int> _syndrome, int _maxNumberOfIterations = 30);
  std::vector<int> decoderBealivePropagation(std::vector<int> _message, std::vector<int> _syndrome, int _maxNumberOfIterations = 30);
  const std::vector<std::vector<int>> &getH() { return H; }
};
}  // namespace LDPC

#endif