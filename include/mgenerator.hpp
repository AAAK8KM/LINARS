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
    
        for (uint32_t i=0;i<n;i++)
            for (uint32_t j=0;j<m;j++)
            {
                M[i*m+j,i*m+j]=4;
                if (0<=i*m+j+1 && i*m+j+1<m) [[likely]]
                    M[i*m+j,i*m+j+1]=-1;
                if (0<=i*m+j-1 && i*m+j-1<m) [[likely]]
                    M[i*m+j,i*m+j-1]=-1;
                if (0<=i*m+j+m && i*m+j+m<m) [[likely]]
                    M[i*m+j,i*m+j+m]=-1;
                if (0<=i*m+j-m && i*m+j-m<m) [[likely]]
                    M[i*m+j,i*m+j-m]=-1;
            }
        return Mtype(M);
    }
    
}


#endif