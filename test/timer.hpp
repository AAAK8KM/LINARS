#ifndef timer_hpp__
#define timer_hpp__

#include <utility>
#include <chrono>

template<typename dret, typename... dargs>
std::pair<dret, double> ms_timer(dret(&func)(dargs...), const dargs&...  args)
{
    auto start = std::chrono::high_resolution_clock::now();
    auto res = func(args...);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()/ 1000.0;
    return std::make_pair(res, duration);
}

#endif