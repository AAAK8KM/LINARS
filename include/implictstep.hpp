#ifndef implictstep_h__
#define implictstep_h__

#include "imatrix.hpp"
#include "t2m.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <chebishev.hpp>
#include <ostream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace LINARS {

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> SSORStep(const Mtype& A, const VMatrix<dtype>& b,const VMatrix<dtype>& prev, dtype w=1)
    {
        if (b.size().first!=prev.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> step(prev);
        Vector<dtype> x(b.size().first);
        std::vector<std::tuple<uint32_t,uint32_t,dtype>> promiseU, promiseL; //promised for whom, required, const
        promiseU.reserve(b.size().first*4); //4 is a funny constant
        promiseL.reserve(b.size().first*4);
        for (auto [i,j,c]: A)
                if (i<j) promiseU.emplace_back(i,j,c);
                else if (i>j) promiseL.emplace_back(i,j,c); //will be fast because of smart alloc
        std::sort(promiseL.begin(),promiseL.end(),std::less<decltype(promiseL.back())>());
        std::sort(promiseU.begin(),promiseU.end(),std::greater<decltype(promiseU.back())>());
        std::array<std::reference_wrapper<decltype(promiseL)>, 2> p={std::ref(promiseL), std::ref(promiseU)};
        for (uint32_t si=0;si<b.size().second;si++)
        {
            for (uint32_t s=0;s<2;s++){
                x=w*b[si];
                for (const auto& [i,j,c]: p[(s+1)&1].get())
                    x[i]-=step[si][j]*c*w;
                for (uint32_t i=0;i<A.size().first;i++)
                {
                    x[i]/=A.ge(i,i);
                    x[i]-=step[si][i]*(w-1);
                }
                for (const auto& [i,j,c]: p[s].get())
                {
                    x[i]-=w*x[j]*c/A.ge(i,i);
                }
                step[si]=x;
            }
        }
        return  step;
    }

}


#endif