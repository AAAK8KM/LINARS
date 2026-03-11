#ifndef implictsolver_h__
#define implictsolver_h__

#include "imatrix.hpp"
#include "t2m.hpp"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace LINARS {

    constexpr uint32_t preset_max_iter=1000;
    constexpr long double preset_max_r=1e-6;

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> JakobiSolver(const Mtype& A, const VMatrix<dtype>& b, uint32_t max_iter=preset_max_iter, dtype max_r=preset_max_r)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();
        uint32_t iter=0;
        Vector<dtype> x(b.size().first);
        Vector<dtype> r;
        while (iter++<max_iter && mes_r>max_r) {
            mes_r=0;
            for (uint32_t si=0;si<b.size().second;si++)
            {
                x=b[si];
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
    VMatrix<dtype> GaussZeidelSolver(const Mtype& A, const VMatrix<dtype>& b, uint32_t max_iter=preset_max_iter, dtype max_r=preset_max_r)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();;
        uint32_t iter=0;
        Vector<dtype> x(b.size().first);
        Vector<dtype> r;
        std::vector<std::tuple<uint32_t,uint32_t,dtype>> promise; //promised for whom, required, const
        promise.reserve(b.size().first*4); //4 is a funny constant
        while (iter++<max_iter && mes_r>max_r) {
            mes_r=0;
            for (uint32_t si=0;si<b.size().second;si++)
            {
                x=b[si];
                for (auto [i,j,c]: A)
                    if (i<j) x[i]-=sol[si][j]*c;
                    else if (i!=j) promise.emplace_back(i,j,c); //will be fast because of smart alloc
                
                std::sort(promise.begin(),promise.end(),std::greater<decltype(promise.back())>());
                
                for (uint32_t i=0;i<A.size().first;i++)
                    x[i]/=A.gev(i,i);

                for (;promise.size()>0;)
                {
                    auto [i,j,c] = promise.back();
                    //std::cout<<i<<" "<<j<<std::endl;
                    x[i]-=x[j]*c/A.gev(i,i);
                    promise.pop_back();
                }
                r=A*x-b;
                //std::cout<<"err="<<(r|r)<<" "<<iter<<std::endl<<x<<std::endl;
                mes_r=std::max(mes_r,r|r);
                sol[si]=x;
            }
        }
        return  sol;
    }


    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> SimpleSolver(const Mtype& A, const VMatrix<dtype>& b, dtype tau, uint32_t max_iter=preset_max_iter, dtype max_r=preset_max_r)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();
        uint32_t iter=0;
        while (iter++<max_iter && mes_r>max_r) {
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