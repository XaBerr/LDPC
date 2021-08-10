#ifndef LDPC_GALLAGER_H
#define LDPC_GALLAGER_H

#include <vector>

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

 public:
  Gallager(int _n, int _k);
  // const std::vector<int> &decode(const std::vector<int> &_codeword);
  // const std::vector<int> &encode(const std::vector<int> &_message);
  const std::vector<int> &getSyndrome(const std::vector<int> &_message);
  // const std::vector<int> &decode(const std::vector<int> &_message, const std::vector<int> &_syndrome);
  std::vector<int> decoderBitFlip(std::vector<int> _message, std::vector<int> _syndrome, int _maxNumberOfIterations = 10);
  const std::vector<std::vector<int>> &getH() { return H; }
};
}  // namespace LDPC

#endif