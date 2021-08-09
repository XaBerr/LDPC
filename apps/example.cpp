#include <LDPC.hpp>
#include <iostream>
#include <random>

using namespace LDPC;
static Uniform uniform;

void channel(std::vector<int>& _message, std::vector<int>& _syndrome, float _pError) {
  for (auto& var : _message)
    if (uniform() < _pError)
      var = (var == 0) ? 1 : 0;

  for (auto& var : _syndrome)
    if (uniform() < _pError)
      var = (var == 0) ? 1 : 0;
}

int main(int argc, char const* argv[]) {
  const int n = 30;
  const int k = 21;
  std::vector<int> messageInput(k);
  std::vector<int> syndromeInput(n - k);
  std::vector<int> messageChannel(k);
  std::vector<int> syndromeChannel(n - k);
  std::vector<int> messageOutput(k);
  std::vector<int> syndromeOutput(n - k);

  for (auto& var : messageInput)
    var = (uniform() < 0.7) ? 0 : 1;

  Gallager ldpc(n, k);
  syndromeInput   = ldpc.getSyndrome(messageInput);
  messageChannel  = messageInput;
  syndromeChannel = syndromeInput;
  channel(messageChannel, syndromeChannel, 0.1);
  messageOutput = ldpc.decode(messageChannel, syndromeChannel);

  std::cout << "Message input: ";
  for (auto var : messageInput)
    std::cout << var;

  std::cout << "\nMessage output: ";
  for (auto var : messageOutput)
    std::cout << var;
}
