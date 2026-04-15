#ifndef timer_hpp__
#define timer_hpp__

#include <utility>
#include <chrono>

template<typename dfunc, typename... dargs>
auto ms_timer(dfunc& func, const dargs&... args)
{
    auto start = std::chrono::system_clock::now();
    auto res = func(args...);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    return std::make_pair(res, duration.count());
}

#endif