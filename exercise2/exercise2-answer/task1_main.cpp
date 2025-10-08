// main.cpp
#include "task1_timer.hpp"
#include <thread>

int main() {

  Timer<> timer;
  std::this_thread::sleep_for(std::chrono::seconds(1));
  return 0;
}
