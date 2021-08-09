#ifndef LDPC_UNIFORM_H
#define LDPC_UNIFORM_H

#include <random>
#include <ctime>

namespace LDPC {
class Uniform {
 private:
  std::mt19937 mt;
  std::uniform_real_distribution<double> dist;

 public:
  Uniform(double lower = 0.0, double upper = 1.0) : mt(std::time(nullptr)), dist(lower, upper) {}
  double operator()() { return dist(mt); }
};
}  // namespace LDPC

#endif