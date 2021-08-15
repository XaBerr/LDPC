#include <LDPC.hpp>
#include <iostream>
#include <random>

using namespace LDPC;
static Uniform uniform;

void channel(std::vector<uint8_t>& _codeword, float _pError) {
  for (auto& var : _codeword)
    if (uniform() < _pError)
      var = (var == 0) ? 1 : 0;
}

void printHexaVector(const std::vector<uint8_t>& _codeword) {
  unsigned int out = 0;
  for (int z = 0; z < _codeword.size(); z++) {
    out |= (_codeword[z] << (z % 4));
    if (z > 0 && (z % 4 == 0)) {
      std::cout << " " << std::uppercase << std::hex << out;
      out = 0;
    }
  }
  if (_codeword.size() % 4)
    std::cout << " " << std::uppercase << std::hex << out;
}

void printBinVector(const std::vector<uint8_t>& _codeword) {
  for (auto var : _codeword)
    std::cout << (int)var;
}

void printMatrix(const std::vector<std::vector<uint8_t>>& _matrix) {
  for (auto row : _matrix) {
    for (auto col : row)
      std::cout << (int)col;
    std::cout << std::endl;
  }
}

int main(int argc, char const* argv[]) {
  const int n = 100;
  const int k = 10;
  const int m = n - k;
  std::vector<uint8_t> messageInput(k);
  std::vector<uint8_t> syndrome(m);
  std::vector<uint8_t> codewordInput(m);
  std::vector<uint8_t> codewordChannel(m);
  std::vector<uint8_t> codewordOutput(m);

  for (auto& var : messageInput)
    var = (uniform() < 0.7) ? 0 : 1;

  Gallager ldpc(n, k, 10);
  codewordInput   = ldpc.getCodeword(messageInput);
  codewordChannel = codewordInput;
  channel(codewordChannel, 0.05);
  codewordOutput = ldpc.decoderBealivePropagation(codewordChannel);

  syndrome = ldpc.checkSyndrome(codewordInput);
  std::cout << "\Syndrome check: ";
  printHexaVector(syndrome);

  // std::cout << "\nH matrix: " << ldpc.H.size() << " " << ldpc.H[5].size() << std::endl;
  // printMatrix(ldpc.H);

  // std::cout << "\nHRowEchelon matrix: " << ldpc.HRowEchelon.size() << " " << ldpc.HRowEchelon[5].size() << std::endl;
  // printMatrix(ldpc.HRowEchelon);

  // std::cout << "\nG matrix: " << ldpc.G.size() << " " << ldpc.G[5].size() << std::endl;
  // printMatrix(ldpc.G);

  // std::cout << "\nMessage: ";
  // printHexaVector(messageInput);

  // std::cout << "\nEncoder: ";
  // printHexaVector(codewordInput);

  // std::cout << "\nChannel: ";
  // printHexaVector(codewordChannel);

  // std::cout << "\nDecoder: ";
  // printHexaVector(codewordOutput);
}
