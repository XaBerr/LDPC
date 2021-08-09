#include <LDPC.hpp>
#include <iostream>
#include <random>

using namespace LDPC;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dist(0.0, 1.0);

void channel(std::vector<int>& _message, std::vector<int>& _syndrome, float _pError) {
  for (auto& var : _message)
    if (dist(gen) < _pError)
      var = (var == 0) ? 1 : 0;

  for (auto& var : _syndrome)
    if (dist(gen) < _pError)
      var = (var == 0) ? 1 : 0;
}

int main(int argc, char const* argv[]) {
  const int n = 10;
  const int k = 7;
  std::vector<int> messageInput(k);
  std::vector<int> syndromeInput(n - k);
  std::vector<int> messageChannel(k);
  std::vector<int> syndromeChannel(n - k);
  std::vector<int> messageOutput(k);
  std::vector<int> syndromeOutput(n - k);

  for (auto& var : messageInput)
    var = (dist(gen) < 0.7) ? 0 : 1;

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
