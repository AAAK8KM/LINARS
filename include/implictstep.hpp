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
    std::array<std::vector<std::tuple<uint32_t,uint32_t,dtype>>, 2> SSORPrep(const Mtype& A)
    {
        std::array<std::vector<std::tuple<uint32_t,uint32_t,dtype>>,2> p; //promised for whom, required, const
        std::vector<std::tuple<uint32_t,uint32_t,dtype>>& promiseL = p[0], &promiseU = p[1];
        promiseU.reserve(A.size().second*4); //4 is a funny constant
        promiseL.reserve(A.size().second*4);
        for (auto [i,j,c]: A)
                if (i<j) {if (c!=0) promiseU.emplace_back(i,j,c);}
                else if (i>j) if (c!=0) promiseL.emplace_back(i,j,c); //will be fast because of smart alloc
        std::sort(promiseL.begin(),promiseL.end(),std::less<>{});
        std::sort(promiseU.begin(),promiseU.end(),std::greater<>{});
        return p;
    }

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> SSORStep(const Mtype& A, const VMatrix<dtype>& b,const VMatrix<dtype>& prev,  const std::array<std::vector<std::tuple<uint32_t,uint32_t,dtype>>, 2>& p, dtype w=1)
    {
        if (b.size().first!=prev.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> step(prev);
        Vector<dtype> x(b.size().first);
        for (uint32_t si=0;si<b.size().second;si++)
        {
            for (uint32_t s=0;s<2;s++){
                x=w*b[si];
                for (const auto& [i,j,c]: p[(s+1)&1])
                    x[i]-=step[si][j]*(c*w);
                /*std::for_each(__pstl::execution::par_unseq, p[(s+1)&1].begin(), p[(s+1)&1].end(), [&x,&step,w,si](const auto& xv){
                    const auto& [i,j,c]=xv;
                   
                });*/
                for (uint32_t i=0;i<A.size().first;i++)
                {
                    x[i]/=A.ge(i,i);
                    x[i]-=step[si][i]*(w-1);
                }
                for (const auto& [i,j,c]: p[s])
                    x[i]-=x[j]*(w*c/A[i,i]);
                step[si]=x;
            }
        }
        return  step;
    }

}


#endif