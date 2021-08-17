#ifndef LDPC_LDPC_H
#define LDPC_LDPC_H

#include <vector>
#include <math.h>
#include <stdint.h>

namespace LDPCm {
class LDPC {
 public:
  const size_t n;
  const size_t k;
  const size_t m;
  const size_t numberOfRowsInH;
  std::vector<std::vector<uint8_t>> H;
  std::vector<std::vector<uint8_t>> HRowEchelon;
  std::vector<std::vector<uint8_t>> G;

 protected:
  std::vector<std::vector<size_t>> codewordsLinks;
  std::vector<std::vector<size_t>> paritiesLinks;
  std::vector<uint8_t> message;
  std::vector<uint8_t> syndrome;
  std::vector<uint8_t> codeword;

  static inline float loglikelihood(float d) { return log(d / (1 - d)); }
  // static inline float atanh(float f) { return log((1 + f) / (1 - f)); }
  static inline float phy(float x) { return log((exp(x) + 1) / (exp(x) - 1)); }

 public:
  LDPC(size_t _n, size_t _k, size_t _numberOfRowsInH);
  const std::vector<std::vector<uint8_t>> &generateH() { return H; }
  const std::vector<std::vector<uint8_t>> &generateHRowEchelon();
  const std::vector<std::vector<uint8_t>> &generateG();
  void generateLinks();
  void generateHG(int _maxNumberOfIterations = 100) {}

  const std::vector<uint8_t> &checkSyndrome(const std::vector<uint8_t> &_codeword);
  const std::vector<uint8_t> &getCodeword(const std::vector<uint8_t> &_message);
  std::vector<uint8_t> decoderBitFlip(std::vector<uint8_t> _message, std::vector<uint8_t> _syndrome, size_t _maxNumberOfIterations = 30);
  std::vector<uint8_t> decoderBealivePropagation(std::vector<uint8_t> _codeword, size_t _maxNumberOfIterations = 30);
};
}  // namespace LDPCm

#endif