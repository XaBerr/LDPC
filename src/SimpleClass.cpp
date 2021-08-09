#include <LDPC/SimpleClass.hpp>

using namespace LDPC;

SimpleClass::SimpleClass()
    : number{0} {}

void SimpleClass::setNumber(int _number) {
  number = _number;
}
int SimpleClass::getNumber() const {
  return number;
}