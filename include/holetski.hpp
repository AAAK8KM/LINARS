#ifndef holetski_hpp__
#define holetski_hpp__

#include "t2m.hpp"
#include "matrix.hpp"
#include "mcsr.hpp"
#include "mdok.hpp"
#include "mvector.hpp"
#include "types.hpp"
#include <cmath>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <utility>

namespace LINARS {

    template<typename dtype>
    Vector<dtype> MulRLLT(const MCSR<dtype> &L, const Vector<dtype> &v)
    {
        Vector<dtype> r(v.size());
        for (uint32_t i=0;i<L.size().first;i++)
        {
            r[i]=v[i];
            for (auto [j,c]: L.gRowSpan(i))
                if (i!=j) r[i]-=r[j]*c;
            r[i]=r[i]/L[i,i];
        }
        //std::cout<<r;
        for (uint32_t i=L.size().first-1;i<L.size().first;i--)
        {
            r[i]=r[i]/L[i,i];
            for (auto [j,c]: L.gRowSpan(i))
                if (i!=j) r[j]-=r[i]*c;
        }
        return r;
    }


    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    MCSR<dtype> holetski(const Mtype& M)
    {
        if (M.size().first!=M.size().second) throw std::invalid_argument("holetski is only for square matrixes.");
        uint32_t n=M.size().first;
        MDOK<dtype> L(M.size());

        // L is L^T
        //.00 ...
        //..0 0..
        //... 00.
        //L[i,j] = M[i,j]-sum^{i-1}_k=1 (L[i,k]*L[j,k]) i<=j
        //L[i,i] = M[i,i]-sum^{i-1}_k=1 (L[i,k]*L[i,k]) for [i,i] we ned all prev horiz 

        //std::cout<<n<<" s\n";

        for (uint32_t i=0;i<n;i++)
        {
            for (uint32_t j=i;j<n;j++)
            {
                if (M.gev(i,j)!=0)
                {
                    L[j,i]=M[i,j];
                    //std::cout<<i<<" "<<j<<std::endl;
                    for (uint32_t k=0;k<i;k++)
                        L[j,i]-=L[i,k]*L[j,k];
                    //std::cout<<L[j,i]<<std::endl;
                    if (i==j) L[j,i]=std::sqrt(L[j,i]);
                    else L[j,i]/=L[i,i];
                    //std::cout<<L[j,i]<<std::endl;
                }
            }
        } 
        std::cout<<"e\n";
        return MCSR<dtype>(L);
    }

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    MCSR<dtype> holetski(const Mtype& M, const Matrix<uint32_t>& LVL, uint32_t p)
    {
        if (M.size().first!=M.size().second) throw std::invalid_argument("holetski is only for square matrixes.");
        if (LVL.size().first!=M.size().second) throw std::invalid_argument("LVL wrong size.");
        if (LVL.size().first!=LVL.size().second) throw std::invalid_argument("LVL is not square.");
        uint32_t n=M.size().first;
        MDOK<dtype> L(M.size());

        for (uint32_t j=0;j<n;j++)
        {
            for (uint32_t i=0;i<=j;i++)
            {
                if (M.gev(i,j)!=0 && LVL[i,j]<=p)
                {
                    L[j,i]=M[i,j];
                    for (uint32_t k=0;k<j-1;k++)
                        L[j,i]-=L[i,k]*L[j,k];
                    if (i==j) L[j,i]=std::sqrt(L[j,i]);
                    else L[j,i]/=L[i,i];
                }
            }
        } 

        return MCSR<dtype>(L);
    }
}


#endif