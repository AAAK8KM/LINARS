#ifndef chebishev_h__
#define chebishev_h__

#include <array>
#include <cassert>
#include <cstdint>
#include <numbers>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace LINARS{

template <uint32_t N>
consteval uint32_t log2()
{
    uint32_t n=N, i=0;
    while (n>1) 
    {
        if (n&1) assert(false);
        n/=2;
        i++;
    }
    return i;
}

template <uint32_t N>
consteval std::array<uint32_t, N> ChebPermutations()
{
    uint32_t j=0;
    uint32_t mask;
    std::array<uint32_t, N> res;
    assert(log2<N>());
    for (uint32_t i=0;i<N;i++)
    {
        res[i]=j;
        mask =  N>>1; 
        while (mask & j) {
            j ^= mask;
            mask >>= 1;
        }
        j ^= mask;
    }
    return res;
}


template <uint32_t N, typename dtype> 
consteval std::array<dtype, N> ChebRoots()
{
    dtype cos0,sin0;
    std::array<dtype, N> ccos, csin;
    ccos[0]=std::cos(std::numbers::pi_v<long double>/(2*N));
    csin[0]=std::sin(std::numbers::pi_v<long double>/(2*N));
    cos0=std::cos(std::numbers::pi_v<long double>/N);
    sin0=std::sin(std::numbers::pi_v<long double>/N);
    for (uint32_t i=1;i<N;i++)
    {
        ccos[i]=ccos[i-1]*cos0-csin[i-1]*sin0;
        csin[i]=csin[i-1]*cos0+ccos[i-1]*sin0;
    }
    return ccos;
}

template <uint32_t N, typename dtype> 
consteval std::array<dtype, N> ChebRootsPermutaton()
{
    std::array<uint32_t, N> p=ChebPermutations<N>();
    std::array<dtype, N> roots=ChebRoots<N, dtype>();
    for (uint32_t i=0;i<N;i++)
        if (i<p[i]) std::swap(roots[i], roots[p[i]]);
    return roots;
}

//transforms from [-1, 1] to [a, b]
template <uint32_t N, typename dtype> 
std::array<dtype, N> LinTrasf(const std::array<dtype, N>& v, dtype a, dtype b)
{
    if (a>b) throw std::runtime_error("wrong parametrs for linar transform");
    std::array<dtype, N> res(v);
    for (auto& c: res)
        c=c*(b-a)/2+(a+b)/2;
    return res;
}

}

#endif