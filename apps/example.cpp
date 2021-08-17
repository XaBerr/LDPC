#include <LDPCm.hpp>
#include <iostream>
#include <random>

using namespace LDPCm;
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
  const int n = 8;
  const int k = 4;
  const int m = n - k;
  std::vector<uint8_t> messageInput(k, 1);
  std::vector<uint8_t> syndrome(m);
  std::vector<uint8_t> codewordInput(m);
  std::vector<uint8_t> codewordChannel(m);
  std::vector<uint8_t> codewordOutput(m);

  Gallager ldpc(12, 3, 3);
  ldpc.H = {
      {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1},

      {1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0},
      {0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1},
      {0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0},

      {1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0},
      {0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0},
      {0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1}};
  // ldpc.generateHG();
  // ldpc.generateH();
  // ldpc.generateHRowEchelon();
  // ldpc.generateG();

  // Gallager ldpc(8, 4, 2);
  // ldpc.H = {
  //     {1, 1, 1, 1, 0, 0, 0, 0},
  //     {0, 0, 0, 0, 1, 1, 1, 1},

  //     {0, 0, 1, 0, 1, 0, 1, 1},
  //     {1, 1, 0, 1, 0, 1, 0, 0}};
  // std::cout << "\nH matrix: " << ldpc.H.size() << " " << ldpc.H[5].size() << std::endl;
  // printMatrix(ldpc.H);
  // ldpc.generateHRowEchelon();
  // ldpc.generateG();

  // for (auto& var : messageInput)
  //   var = (uniform() < 0.7) ? 0 : 1;

  codewordInput = ldpc.getCodeword(messageInput);
  syndrome      = ldpc.checkSyndrome(codewordInput);
  // codewordChannel = codewordInput;
  // channel(codewordChannel, 0.05);
  // codewordOutput = ldpc.decoderBealivePropagation(codewordChannel);

  // std::cout << "\Syndrome check: ";
  // printHexaVector(syndrome);

  std::cout << "\nH matrix: " << ldpc.H.size() << " " << ldpc.H[0].size() << std::endl;
  printMatrix(ldpc.H);

  std::cout << "\nHRowEchelon matrix: " << ldpc.HRowEchelon.size() << " " << ldpc.HRowEchelon[0].size() << std::endl;
  printMatrix(ldpc.HRowEchelon);

  // std::cout << "\nG matrix: " << ldpc.G.size() << " " << ldpc.G[0].size() << std::endl;
  // printMatrix(ldpc.G);

  // std::cout << "\nMessage: ";
  // printBinVector(messageInput);

  // std::cout << "\nEncoder: ";
  // printBinVector(codewordInput);

  // std::cout << "\nSyndrome: ";
  // printBinVector(syndrome);

  // std::cout << "\nChannel: ";
  // printHexaVector(codewordChannel);

  // std::cout << "\nDecoder: ";
  // printHexaVector(codewordOutput);
}
