#ifndef LDPC_SIMPLECLASS_H
#define LDPC_SIMPLECLASS_H

namespace LDPCm {
class SimpleClass {
  int number;

 public:
  SimpleClass();
  void setNumber(int _number);
  int getNumber() const;
};
}  // namespace LDPCm

#endif