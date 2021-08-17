#ifndef LDPC_UNIFORM_H
#define LDPC_UNIFORM_H

#include <random>
#include <ctime>

namespace LDPCm {
class Uniform {
  std::mt19937 generator;
  std::uniform_real_distribution<double> distribution;

 public:
  Uniform(double _lower = 0.0, double _upper = 1.0) : generator(std::time(nullptr)), distribution(_lower, _upper) {}
  std::mt19937& getGenerator() { return generator; }
  double operator()() { return distribution(generator); }
};
}  // namespace LDPCm

#endif