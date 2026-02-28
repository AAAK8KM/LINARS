#ifndef tomas_hpp__
#define tomas_hpp__

#include "matrix.hpp"
#include <cstdint>
#include <stdexcept>

namespace LINARS {

template<typename dtype>
class TomasSolver
{
    public:
    static Matrix<dtype> Solve(const M3diag<dtype>& A,const  Matrix<dtype>& d)
    {
        if (A.size().first!=d.size().first) throw std::runtime_error("Invalid system of equations");
        std::uint32_t n=A.size().first,m=d.size().second;
        Matrix<dtype> res(n,m), p(n,m), q(n,m);
        #define a(i) A.ge(i,i-1)
        #define b(i) A.ge(i,i)
        #define c(i) A.ge(i,i+1)
        for (std::uint32_t i=0;i<m;i++)
        {
            q.ge(1,i)=d.ge(0,i)/b(0);
            p.ge(1,i)=-c(0)/b(0);
        }
        for (std::uint32_t i=1;i<n-1;i++)
            for (std::uint32_t j=0;j<m;j++)
            {
                if (i<n-1)
                    p.ge(i+1,j)=-c(i)                   /(b(i)+a(i)*p.ge(i,j));
                q.ge(i+1,j)=(d.ge(i,j)-a(i)*q.ge(i,j))/(b(i)+a(i)*p.ge(i,j));
            }
        for (std::uint32_t i=0;i<m;i++)
            res.ge(n-1,i)=(d.ge(n-1,i)-a(n-1)*q.ge(n-1,i))/(b(n-1)+a(n-1)*p.ge(n-1,i));
        #undef a
        #undef b
        #undef c
        //std::cout<<p<<q;
        for (std::uint32_t i=n-2;i>0;i--)
            for (std::uint32_t j=0;j<m;j++)
                res.ge(i,j)=q.ge(i+1,j)+res.ge(i+1,j)*p.ge(i+1,j);
        for (std::uint32_t j=0;j<m;j++)
            res.ge(0,j)=q.ge(1,j)+res.ge(1,j)*p.ge(1,j);
        return res;
    }
};


#ifndef Dtomas
extern template class TomasSolver<double>;
extern template class TomasSolver<float>;
#endif

}

#endif