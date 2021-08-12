#ifndef LDPC_GALLAGER_H
#define LDPC_GALLAGER_H

#include <vector>
#include <math.h>
#include <stdint.h>

namespace LDPC {
class Gallager {
  const int n;
  const int k;
  const int m;
  std::vector<std::vector<uint8_t>> H;
  std::vector<std::vector<uint8_t>> codewordsLinks;
  std::vector<std::vector<uint8_t>> paritiesLinks;
  std::vector<uint8_t> message;
  std::vector<uint8_t> syndrome;
  std::vector<uint8_t> codeword;

  int validateN(int _n, int _k);
  int validateK(int _n, int _k);
  static inline float loglikelihood(float d) { return log(d / (1 - d)); }
  // static inline float atanh(float f) { return log((1 + f) / (1 - f)); }
  static inline float phy(float x) { return log((exp(x) + 1) / (exp(x) - 1)); }

 public:
  Gallager(int _n, int _k);
  // const std::vector<uint8_t> &decode(const std::vector<uint8_t> &_codeword);
  // const std::vector<uint8_t> &encode(const std::vector<uint8_t> &_message);
  const std::vector<uint8_t> &getSyndrome(const std::vector<uint8_t> &_message);
  // const std::vector<uint8_t> &decode(const std::vector<uint8_t> &_message, const std::vector<uint8_t> &_syndrome);
  std::vector<uint8_t> decoderBitFlip(std::vector<uint8_t> _message, std::vector<uint8_t> _syndrome, int _maxNumberOfIterations = 30);
  std::vector<uint8_t> decoderBealivePropagation(std::vector<uint8_t> _message, std::vector<uint8_t> _syndrome, int _maxNumberOfIterations = 30);
  const std::vector<std::vector<uint8_t>> &getH() { return H; }
  std::vector<std::vector<uint8_t>> eliminationGaussJordan();
  std::vector<std::vector<uint8_t>> gaussReduce(std::vector<std::vector<uint8_t>> H);
};
}  // namespace LDPC

#endif