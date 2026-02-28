#ifndef t2m_h__
#define t2m_h__

#include <cstddef>
#include <istream>
#include "imatrix.hpp"
#include "matrix.hpp"


template<typename dtype>
std::istream& operator>>(std::istream& in, LINARS::Matrix<dtype>& m)
{
    auto p=m.size();
    for (std::size_t i=0;i<p.first;i++)
        for (std::size_t j=0;j<p.second;j++)
            in>>m.ge(i, j);
    return  in;
}

template<typename dtype>
std::ostream& operator<<(std::ostream& o,const LINARS::IMatrix<dtype>& m)
{
    auto p=m.size();
    for (std::size_t i=0;i<p.first;i++)
    {
        for (std::size_t j=0;j<p.second;j++)
            o<<m.gev(i, j)<<" ";
        o<<std::endl;
    }
    return  o;
}

#endif