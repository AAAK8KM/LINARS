#ifndef timer_hpp__
#define timer_hpp__

#include <cstdint>
#include <utility>
#include <chrono>

template<uint32_t N=1,typename dfunc, typename... dargs>
auto ms_timer(dfunc& func, const dargs&... args)
{
    auto start = std::chrono::system_clock::now();
    auto res = func(args...);
    for (uint32_t i=1;i<N;i++)
        res = func(args...);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    return std::make_pair(res, duration.count()/N);
}

#endif