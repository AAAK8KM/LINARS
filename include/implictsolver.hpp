#ifndef implictsolver_h__
#define implictsolver_h__

#include "imatrix.hpp"
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace LINARS {

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> YakobiSolver(const Mtype& A, const VMatrix<dtype>& b, uint32_t max_iter=1000, dtype max_r=1e-6)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();
        uint32_t iter=0;
        Vector<dtype> x(b.size().first);
        Vector<dtype> r;
        while (iter++<max_iter && mes_r<max_r) {
            mes_r=0;
            for (uint32_t si=0;si<b.size().second;si++)
            {
                x=b;
                for (auto [i,j,c]: A)
                    if (i!=j) x[i]-=sol[si][j]*c;
                for (uint32_t i=0;i<A.size().first;i++)
                    x[i]/=A.gev(i,i);
                r=A*x-b;
                mes_r=std::max(mes_r,r|r);
                sol[si]=x;
            }
        }
        return  sol;
    }

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> GaussZeidelSolver(const Mtype& A, const VMatrix<dtype>& b, uint32_t max_iter=1000, dtype max_r=1e-6)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();
        uint32_t iter=0;
        Vector<dtype> x(b.size().first);
        Vector<dtype> r;
        std::vector<std::tuple<uint32_t,uint32_t,dtype>> promise; //promised for whom, required, const
        promise.reserve(b.size().first*4); //4 is a funny constant
        while (iter++<max_iter && mes_r<max_r) {
            mes_r=0;
            for (uint32_t si=0;si<b.size().second;si++)
            {
                x=b;
                for (auto [i,j,c]: A)
                    if (i<j) x[i]-=sol[si][j]*c;
                    else if (i!=j) promise.emplace_back(i,j,c); //will be fast because of smart alloc
                
                std::sort(promise.begin(),promise.end());
                
                for (;promise.size()>0;)
                {
                    auto [i,j,c] = promise.pop_back();
                    x[i]-=c*x[j];
                }

                for (uint32_t i=0;i<A.size().first;i++)
                    x[i]/=A.gev(i,i);
                r=A*x-b;
                mes_r=std::max(mes_r,r|r);
                sol[si]=x;
            }
        }
        return  sol;
    }


    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> SimpleSolver(const Mtype& A, const VMatrix<dtype>& b, dtype tau, uint32_t max_iter=1000, dtype max_r=1e-6)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();
        uint32_t iter=0;
        while (iter++<max_iter && mes_r<max_r) {
            mes_r=0;
            for (uint32_t i=0;i<b.size().second;i++)
            {
                Vector<dtype> r=A*sol[i]-b;
                sol[i]=sol[i]-r*tau;
                mes_r=std::max(mes_r,r|r);
            }
        }
        return sol;
    }


}


#endif