#pragma once
#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <iostream>
#include <string_view>

// Timer class for measuring elapsed time
template <typename T = std::chrono::milliseconds>
struct Timer {
    using Clock = std::chrono::high_resolution_clock;
    using Duration = T;

    Clock::time_point start_time_wall;

    Timer() {
        start_time_wall = Clock::now(); // Start timing
    }

    ~Timer() {
        // Print elapsed time
        std::cout << *this << std::endl;
    }

    void reset() {
        start_time_wall = Clock::now(); // Reset start time
    }

    T elapsed() const {
        // Return elapsed time
        return std::chrono::duration_cast<Duration>(Clock::now() - start_time_wall);
    }

    friend std::ostream &operator<<(std::ostream &os, const Timer &timer) {
        // Output elapsed time
        os << timer.elapsed().count() << " " << time_unit_name<Duration>();
        return os; // Ensure return value
    }

    template <typename U>
    static constexpr std::string_view time_unit_name() {
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
        else
            return "?";
    }
};

#endif // TIMER_HPP