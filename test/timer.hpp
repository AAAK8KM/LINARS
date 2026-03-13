#ifndef timer_hpp__
#define timer_hpp__

#include <utility>
#include <chrono>

template<typename dfunc, typename... dargs>
auto ms_timer(dfunc& func, const dargs&... args)
{
    auto start = std::chrono::high_resolution_clock::now();
    auto res = func(args...);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()/ 1000.0;
    return std::make_pair(res, duration);
}

#endif