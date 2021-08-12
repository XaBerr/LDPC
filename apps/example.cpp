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
  const int n = 30 * 2;
  const int k = 21 * 2;
  const int m = n - k;
  std::vector<uint8_t> messageInput(k);
  std::vector<uint8_t> syndromeInput(m);
  std::vector<uint8_t> messageChannel(k);
  std::vector<uint8_t> syndromeChannel(m);
  std::vector<uint8_t> messageOutput(k);
  std::vector<uint8_t> syndromeOutput(m);

  for (auto& var : messageInput)
    var = (uniform() < 0.7) ? 0 : 1;

  Gallager ldpc(n, k);
  syndromeInput   = ldpc.getSyndrome(messageInput);
  messageChannel  = messageInput;
  syndromeChannel = syndromeInput;
  channel(messageChannel, syndromeChannel, 0.05);
  messageOutput  = ldpc.decoderBealivePropagation(messageChannel, syndromeChannel);
  syndromeOutput = ldpc.getSyndrome(messageOutput);

  auto H = ldpc.getH();
  std::cout << "\nH matrix: " << H.size() << " " << H[5].size() << std::endl;
  for (auto row : H) {
    for (auto col : row)
      std::cout << (int)col;
    std::cout << std::endl;
  }

  std::cout << "\nEncoder: ";
  for (auto var : messageInput)
    std::cout << (int)var;
  std::cout << " | ";
  for (auto var : syndromeInput)
    std::cout << (int)var;

  std::cout << "\nChannel: ";
  for (auto var : messageChannel)
    std::cout << (int)var;
  std::cout << " | ";
  for (auto var : syndromeChannel)
    std::cout << (int)var;

  std::cout << "\nDecoder: ";
  for (auto var : messageOutput)
    std::cout << (int)var;
  std::cout << " | ";
  for (auto var : syndromeOutput)
    std::cout << (int)var;
}
