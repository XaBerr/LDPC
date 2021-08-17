#ifndef LDPC_GALLAGER_H
#define LDPC_GALLAGER_H

#include <LDPCm/LDPC.hpp>

namespace LDPCm {
class Gallager : public LDPC {
 public:
  Gallager(size_t _n = 16, size_t _k = 5, size_t _weightColumns = 3);
  const std::vector<std::vector<uint8_t>> &generateH();
  void generateHG(int _maxNumberOfIterations = 100);
};
}  // namespace LDPCm

#endif