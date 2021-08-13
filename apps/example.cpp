#include <LDPC.hpp>
#include <iostream>
#include <random>

using namespace LDPC;
static Uniform uniform;

void channel(std::vector<uint8_t>& _message, std::vector<uint8_t>& _syndrome, float _pError) {
  for (auto& var : _message)
    if (uniform() < _pError)
      var = (var == 0) ? 1 : 0;

  for (auto& var : _syndrome)
    if (uniform() < _pError)
      var = (var == 0) ? 1 : 0;
}

int main(int argc, char const* argv[]) {
  const int n = 100;
  const int k = 10;
  const int m = n - k;
  std::vector<uint8_t> messageInput(k);
  std::vector<uint8_t> codewordInput(m);
  std::vector<uint8_t> messageChannel(k);
  std::vector<uint8_t> codewordChannel(m);
  std::vector<uint8_t> messageOutput(k);
  std::vector<uint8_t> codewordOutput(m);

  for (auto& var : messageInput)
    var = (uniform() < 0.7) ? 0 : 1;

  Gallager ldpc(n, k, 10);
  codewordInput   = ldpc.getCodeword(messageInput);
  messageChannel  = messageInput;
  codewordChannel = codewordInput;
  channel(messageChannel, codewordChannel, 0.05);
  codewordOutput = ldpc.decoderBealivePropagation(messageChannel, codewordChannel);

  std::cout << "\nH matrix: " << ldpc.H.size() << " " << ldpc.H[5].size() << std::endl;
  for (auto row : ldpc.H) {
    for (auto col : row)
      std::cout << (int)col;
    std::cout << std::endl;
  }

  std::cout << "\nHRowEchelon matrix: " << ldpc.HRowEchelon.size() << " " << ldpc.HRowEchelon[5].size() << std::endl;
  for (auto row : ldpc.HRowEchelon) {
    for (auto col : row)
      std::cout << (int)col;
    std::cout << std::endl;
  }

  std::cout << "\nG matrix: " << ldpc.G.size() << " " << ldpc.G[5].size() << std::endl;
  for (auto row : ldpc.G) {
    for (auto col : row)
      std::cout << (int)col;
    std::cout << std::endl;
  }

  std::cout << "\nMessage: ";
  for (auto var : messageInput)
    std::cout << (int)var;
  std::cout << "\nEncoder: ";
  for (auto var : codewordInput)
    std::cout << (int)var;

  std::cout << "\nChannel: ";
  for (auto var : codewordChannel)
    std::cout << (int)var;

  std::cout << "\nDecoder: ";
  for (auto var : codewordOutput)
    std::cout << (int)var;
}
