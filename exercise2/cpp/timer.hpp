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

  Timer() { 
    start_time_wall = Clock::now();  // 构造时开始计时
  }
  
  ~Timer() {
    std::cout << *this;  // 析构时自动输出
  }
  
  void reset() { 
    start_time_wall = Clock::now();  // 重置开始时间
  }

  T elapsed() const {
    auto now = Clock::now();
    return std::chrono::duration_cast<T>(now - start_time_wall);
  }

  friend std::ostream &operator<<(std::ostream &os, const Timer &timer) {
    auto duration = timer.elapsed();
    os << "Wall: " << duration.count() << time_unit_name<T>();
    return os;
  }

  template <typename U> static constexpr std::string_view time_unit_name() {
    if constexpr (std::is_same_v<U, std::chrono::hours>)
      return "h";
    else if constexpr (std::is_same_v<U, std::chrono::minutes>)
      return "min";
    else if constexpr (std::is_same_v<U, std::chrono::seconds>)
      return "s";
    else if constexpr (std::is_same_v<U, std::chrono::milliseconds>)
      return "ms";
    else if constexpr (std::is_same_v<U, std::chrono::microseconds>)
      return "us";
    else if constexpr (std::is_same_v<U, std::chrono::nanoseconds>)
      return "ns";
    return "?";
  }
}; // Timer
#endif // TIMER_HPP
