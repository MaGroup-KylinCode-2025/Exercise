#pragma once
#ifndef TIMER_HPP
#define TIMER_HPP
#include <chrono>
#include <iostream>

// Timer class for measuring elapsed time
template <typename T = std::chrono::milliseconds> struct Timer {
  using Clock = std::chrono::high_resolution_clock;
  using Duration = T;

  Clock::time_point start_time_wall;

  Timer() { start_time_wall = Clock::now(); }
  ~Timer() { std::cout << *this << std::endl; }
  void reset() { start_time_wall = Clock::now(); }

  T elapsed() const {
    return std::chrono::duration_cast<Duration>(Clock::now() - start_time_wall);
  }

  friend std::ostream &operator<<(std::ostream &os, const Timer &timer) {
    os << std::fixed << std::setprecision(3) << timer.elapsed().count() << " "
       << Timer::time_unit_name<Duration>();
    return os;
  }

  template <typename U> static constexpr std::string_view time_unit_name() {
    if constexpr (std::is_same_v<U, std::chrono::seconds>)
      return "s";
    else if constexpr (std::is_same_v<U, std::chrono::milliseconds>)
      return "ms";
    else if constexpr (std::is_same_v<U, std::chrono::microseconds>)
      return "us";
    else if constexpr (std::is_same_v<U, std::chrono::nanoseconds>)
      return "ns";
    else if constexpr (std::is_same_v<U, std::chrono::minutes>)
      return "min";
    else if constexpr (std::is_same_v<U, std::chrono::hours>)
      return "h";
    else
      return "?";
  }
}; // Timer
#endif // TIMER_HPP