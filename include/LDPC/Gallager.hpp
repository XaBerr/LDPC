#ifndef LDPC_GALLAGER_H
#define LDPC_GALLAGER_H

#include <vector>
#include <math.h>
#include <stdint.h>

namespace LDPC {
class Gallager {
 public:
  const size_t n;
  const size_t k;
  const size_t m;
  const size_t weightColumns;
  const size_t weightRows;
  std::vector<std::vector<uint8_t>> H;
  std::vector<std::vector<uint8_t>> HRowEchelon;
  std::vector<std::vector<uint8_t>> G;

 private:
  std::vector<std::vector<size_t>> codewordsLinks;
  std::vector<std::vector<size_t>> paritiesLinks;
  std::vector<uint8_t> message;
  std::vector<uint8_t> syndrome;
  std::vector<uint8_t> codeword;

  size_t validateN(size_t _n);
  size_t validateK(size_t _n, size_t _k);
  size_t validateM(size_t _n, size_t _k, size_t _weightColumns);
  static inline float loglikelihood(float d) { return log(d / (1 - d)); }
  // static inline float atanh(float f) { return log((1 + f) / (1 - f)); }
  static inline float phy(float x) { return log((exp(x) + 1) / (exp(x) - 1)); }

 public:
  Gallager(size_t _n = 16, size_t _k = 5, size_t _weightColumns = 3);
  // const std::vector<uint8_t> &decode(const std::vector<uint8_t> &_codeword);
  // const std::vector<uint8_t> &encode(const std::vector<uint8_t> &_message);
  const std::vector<uint8_t> &getSyndrome(const std::vector<uint8_t> &_message);
  // const std::vector<uint8_t> &decode(const std::vector<uint8_t> &_message, const std::vector<uint8_t> &_syndrome);
  std::vector<uint8_t> decoderBitFlip(std::vector<uint8_t> _message, std::vector<uint8_t> _syndrome, size_t _maxNumberOfIterations = 30);
  std::vector<uint8_t> decoderBealivePropagation(std::vector<uint8_t> _message, std::vector<uint8_t> _syndrome, size_t _maxNumberOfIterations = 30);
  const std::vector<std::vector<uint8_t>> &generateH();
  const std::vector<std::vector<uint8_t>> &generateHRowEchelon();
  const std::vector<std::vector<uint8_t>> &generateG();
};
}  // namespace LDPC

#endif