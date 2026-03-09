#ifndef implictsolver_h__
#define implictsolver_h__

#include "imatrix.hpp"
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace LINARS {

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> YakobiSolver(const Mtype& A, const VMatrix<dtype>& b, uint32_t max_iter=1000, dtype max_r=1e-6)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();
        uint32_t iter=0;
        while (iter++<max_iter && mes_r<max_r) {
            mes_r=0;
            for (uint32_t si=0;si<b.size().second;si++)
            {
                Vector<dtype> r(b.size().first);
                for (auto [i,j,c]: A)
                    if (i!=j) r[i]+=sol[si][j]*c;
                for (uint32_t i=0;i<A.size().first;i++)
                    r[i]/=A.gev(i,i);
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