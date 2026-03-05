#ifndef gauss_h__
#define gauss_h__

#include "imatrix.hpp"
#include <cstdint>
#include <stdexcept>
#include <utility>

namespace LINARS {

template<typename dtype>
VMatrix<dtype> ReverseGauss(const IMatrix<dtype>& A, const IMatrix<dtype>& b)
{
    if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
    VMatrix<dtype> res(std::make_pair(A.size().first, b.size().second));
    for (uint32_t row=A.size().first-1;row<A.size().first;row--)
        for (uint32_t sol=0;sol<b.size().second;sol++)
        {
            dtype acc=0;
            for (uint32_t j=A.size().second-1;j>row;j--)
                acc+=A.gev(row, j)*res.gev(j, sol);
            res.ge(row, sol)=(b.gev(row, sol)-acc)/A.gev(row, row);
        }
    return res;
}

template<typename dtype>
Vector<dtype> ReverseGauss(const IMatrix<dtype>& A, const Vector<dtype>& b)
{
    if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
    Vector<dtype> res(A.size().first);
    for (uint32_t row=A.size().first-1;row<A.size().first;row--)
        {
            dtype acc=0;
            for (uint32_t j=A.size().second-1;j>row;j--)
                acc+=A.gev(row, j)*res[j];
            res[row]=(b[row]-acc)/A.gev(row, row);
        }
    return res;
}

}


#endif