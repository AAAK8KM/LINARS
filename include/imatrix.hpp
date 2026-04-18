#ifndef imatrix_hpp__
#define imatrix_hpp__

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>
#include <type_traits>
#include "types.hpp"

namespace LINARS {

template<typename dtype>
class IIterator{
    protected:
        const IMatrix<dtype>* wptr;
    public:
        IIterator(const IMatrix<dtype>& M):wptr(M.get_ptr()){}
        virtual std::tuple<uint32_t,uint32_t,dtype> operator*() = 0;
        virtual void operator++()=0;
        //virtual bool operator!=(const IIterator& it) const = 0;  
}; 

template<typename dtype,typename T>
concept IsMatrix = std::is_base_of_v<IMatrix<dtype>, T>;



template<typename dtype>
class IMatrix
{
    protected:
    const IMatrix<dtype>* get_ptr() const{
        return this;
    }
    public:
        virtual dtype& ge(const uint32_t i, const uint32_t j) = 0;
        virtual const dtype& ge(const uint32_t i, const uint32_t j) const = 0;

        dtype& operator[](const uint32_t i, const uint32_t j)
        {
            return this->ge(i, j);
        };

        const dtype&  operator[](const uint32_t i, const uint32_t j) const
        {
            return this->ge(i, j);
        };

        virtual dtype gev(const uint32_t i, const uint32_t j) const = 0;

        inline virtual const std::pair<uint32_t, uint32_t> size() const = 0;

        template<typename  Mtype2>
        requires IsMatrix<dtype, Mtype2>
        Matrix<dtype> operator*(const Mtype2& B)
        {
            if (this->size().second!=B.size().first) throw std::runtime_error("Matrix and vector has worng sizes. Can not multiply!");
            Matrix<dtype> C(this->size().first,B.size().second);
            for (uint32_t i=0;i<this->size().first;i++)
                for (auto [k,j,c]: B)
                    C[i,j]+=(*this).gev(i,k)*c;
            return C;
        }

        bool operator==(const IMatrix& B) const
        {
            if (this->size().second!=B.size().second) return false;
            if (this->size().first!=B.size().first) return false;

            for (uint32_t i=0;i<this->size().first;i++)
                for (uint32_t j=0;j<this->size().second;j++)
                    if (std::abs(this->gev(i,j)-B.gev(i, j))>1e-9)
                        return false;
            return true;
        }

        bool operator!=(const IMatrix& B) const
        {
            return !(*this==B);
        }

        template<typename T>
        friend class IIterator;

        virtual ~IMatrix(){}
};


template<typename dtype,typename  Mtype>
requires IsMatrix<dtype, Mtype>
Mtype operator*(const Mtype& A,const dtype& value)
{
    Mtype C(A);
    for (auto [i,j,c]: A)
        C.ge(i,j)=c*value;
    return C;
}

template<typename dtype,typename  Mtype>
requires IsMatrix<dtype, Mtype>
Mtype operator/(const Mtype& A,const dtype& value)
{
    Mtype C(A);
    for (auto [i,j,c]: A)
        C.ge(i,j)=c/value;
    return C;
}

template<typename dtype,typename  Mtype>
requires IsMatrix<dtype, Mtype>
Mtype operator*(const dtype& value,const Mtype& A)
{
    return A*value;
}


}

#endif