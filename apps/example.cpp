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

void printVectorHex(const std::vector<uint8_t>& _codeword) {
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

void printVectorBin(const std::vector<uint8_t>& _codeword) {
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
  const int n  = 30;
  const int k  = 22;
  const int m  = n - k;
  const int l  = 8;
  const int wc = 5;
  LDPC ldpc(n, k, l);
  ldpc.H = {{1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0},
            {0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1},
            {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0},
            {0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},
            {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0},
            {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0},
            {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0},
            {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0}};
  // ldpc.H = {
  //     {1, 1, 0},
  //     {0, 1, 1},
  //     {1, 0, 1}};
  // ldpc.H = {
  //     {1, 1, 0, 1, 0, 0},
  //     {0, 1, 1, 0, 1, 0},
  //     {1, 0, 0, 0, 1, 1},
  //     {0, 0, 1, 1, 0, 1}};
  // ldpc.H = {
  //     {1, 1, 0, 1, 1, 0, 0, 1, 0, 0},
  //     {0, 1, 1, 0, 1, 1, 1, 0, 0, 0},
  //     {0, 0, 0, 1, 0, 0, 0, 1, 1, 1},
  //     {1, 1, 0, 0, 0, 1, 1, 0, 1, 0},
  //     {0, 0, 1, 0, 0, 1, 0, 1, 0, 1}};
  ldpc.generateHRowEchelon();
  ldpc.generateG();
  // ldpc.generateHG();

  std::vector<uint8_t> messageInput(k, 1);
  std::vector<uint8_t> syndrome(m);
  std::vector<uint8_t> codewordInput(n);
  std::vector<uint8_t> codewordChannel(n);
  std::vector<uint8_t> codewordOutput(n);

  // std::cout << "\nH matrix: " << ldpc.H.size() << " " << ldpc.H[0].size() << std::endl;
  // printMatrix(ldpc.H);

  // for (auto& var : messageInput)
  //   var = (uniform() < 0.7) ? 0 : 1;

  messageInput       = {0, 0, 1};
  codewordInput      = ldpc.getCodeword(messageInput);
  codewordChannel    = codewordInput;
  codewordChannel[1] = codewordChannel[1] ? 0 : 1;
  codewordChannel[3] = codewordChannel[3] ? 0 : 1;
  // channel(codewordChannel, 0.05);
  codewordOutput = ldpc.decoderBealivePropagation(codewordChannel);

  // std::cout << "\Syndrome check: ";
  // printVectorHex(syndrome);

  std::cout << "\nH matrix: " << ldpc.H.size() << " " << ldpc.H[0].size() << std::endl;
  printMatrix(ldpc.H);

  std::cout << "\nHRowEchelon matrix: " << ldpc.HRowEchelon.size() << " " << ldpc.HRowEchelon[0].size() << std::endl;
  printMatrix(ldpc.HRowEchelon);

  std::cout << "\nG matrix: " << ldpc.G.size() << " " << ldpc.G[0].size() << std::endl;
  printMatrix(ldpc.G);

  std::cout << "\nMessage: ";
  printVectorBin(messageInput);

  std::cout << "\nEncoder: ";
  printVectorBin(codewordInput);
  syndrome = ldpc.checkSyndrome(codewordInput);
  std::cout << " (syndrome) ";
  printVectorBin(syndrome);

  std::cout << "\nChannel: ";
  printVectorBin(codewordChannel);
  syndrome = ldpc.checkSyndrome(codewordChannel);
  std::cout << " (syndrome) ";
  printVectorBin(syndrome);

  std::cout << "\nDecoder: ";
  printVectorBin(codewordOutput);
  syndrome = ldpc.checkSyndrome(codewordOutput);
  std::cout << " (syndrome) ";
  printVectorBin(syndrome);
}
