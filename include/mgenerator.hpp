#ifndef mgenerator_h__
#define mgenerator_h__


#include "imatrix.hpp"
#include "mdok.hpp"
#include <cstdint>


namespace LINARS {
    
    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    Mtype PuassonTask0(uint32_t n, uint32_t m)
    {
        MDOK<dtype> M(m*n,m*n);
        for (uint32_t i=0;i<m;i++)
        {
            M[i,i]=1;
            M[m*n-1-i,m*n-1-i]=1;
        }
        for (uint32_t i=0;i<n;i++)
        {
            M[m*i,m*i]=1;
            M[m*i+m-1,m*i+m-1]=1;
        }
    
        for (uint32_t i=1;i<n-1;i++)
            for (uint32_t j=1;j<m-1;j++)
            {
                M[i*m+j,i*m+j]=4;
                M[i*m+j,i*m+j+1]=-1;
                M[i*m+j,i*m+j-1]=-1;
                M[i*m+j,i*m+j+m]=-1;
                M[i*m+j,i*m+j-m]=-1;
            }
        return Mtype(M);
    }
    
}


#endif