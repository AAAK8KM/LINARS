#ifndef lvlmatrix_hpp__
#define lvlmatrix_hpp__

#include "matrix.hpp"
#include "types.hpp"
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <utility>

namespace LINARS {

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    Matrix<uint32_t> lvlmatrix(const Mtype& M)
    {
        if (M.size().first!=M.size().second) throw std::invalid_argument("lvlmatrix is only for square matrixes.");
        Matrix<uint32_t> LVL(M.size());
        uint32_t n=M.size().first;
        for (auto [i,j,c]: LVL)
            if (M.gev(i,j)==0) LVL[i,j]=n-1;
            else LVL[i,j]=0;
        
        for (uint32_t k=0;k<n;k++)
            for (uint32_t i=k;i<n;i++)
                for (uint32_t j=0;j<n;j++)
                    LVL[i,j]=std::min(LVL[i,k]+LVL[k,j]+1,LVL[i,j]);
        return LVL;
    }

}


#endif