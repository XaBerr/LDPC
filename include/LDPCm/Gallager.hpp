#ifndef LDPC_GALLAGER_H
#define LDPC_GALLAGER_H

#include <LDPCm/LDPC.hpp>

namespace LDPCm {
class Gallager : public LDPC {
 protected:
  size_t validateN(size_t _n);
  size_t validateK(size_t _n, size_t _k);
  size_t validateWeightColumns(size_t _n, size_t _k, size_t _numberOfRowsInH, size_t _weightColumns);

 public:
  const size_t weightColumns;
  const size_t weightRows;
  Gallager(size_t _n, size_t _k, size_t _numberOfRowsInH, size_t _weightColumns);
  const std::vector<std::vector<uint8_t>> &generateH();
  void generateHG(int _maxNumberOfIterations = 100);
};
}  // namespace LDPCm

#endif