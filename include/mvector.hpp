#ifndef vector_hpp__
#define vector_hpp__

#include "types.hpp"
#include "imatrix.hpp"
#include <vector>
#include <cmath>
#include <cstdint>
#include <span>
#include<type_traits>
#include <stdexcept>

namespace LINARS {

template<typename dtype,typename  Mtype>
requires IsMatrix<dtype, Mtype>
Vector<dtype> operator*(const Mtype& A,const Vector<dtype>& B)
{
    if (A.size().second!=B.size().first) throw std::runtime_error("Matrix and vector has worng sizes. Can not multiply!");
    Vector<dtype> C(A.size().first);
    for (auto [i,j,c]: A)
        C[i]+=c*B[j];
    return C;
}

template<typename dtype>
class Vector: public IMatrix<dtype>
{
    private:
        bool transp;
        std::span<dtype> data;
        bool owns;
    public:
        Vector(){}
        Vector(uint32_t n_, bool t=0):transp(t),data(new dtype[n_](0),n_),owns(1){}

        Vector(std::pair<uint32_t, uint32_t> s):transp(s.first<s.second),
        data(new dtype[s.first<s.second?s.second:s.first](0),s.first<s.second?s.second:s.first),owns(1)
        {if ((s.first>s.second?s.second:s.first)!=1) throw std::runtime_error("Wrong vector size to create!");}

        //Vector(std::vector<dtype>& v, bool t=0):transp(t),data(v),owns(1){}

        //Vector(decltype(std::vector<dtype>().begin()) begin, uint32_t size, bool t=0):transp(t),data(begin,size),owns(0){}

        
        Vector(dtype* begin, uint32_t size, bool t=0):transp(t),data(begin,size),owns(0){}

        Vector( uint32_t size, dtype* source, bool t=0):transp(t),data(new dtype[size](0),size),owns(1){
            std::copy(source, source+size, data.begin());
        }

        //                                                                      | я не вижу другого выхода, кроме как
        //                                                                      | полностью переписать логику типов,
        //                                                                      V но пока нет времени заниматься этим
        //Vector(const dtype* begin, uint32_t size, bool t=0):transp(t),data(const_cast<dtype*>(begin),size),owns(0){}
        
        Vector& operator=(const Vector& rhs) 
        {
            if (this->size().first!=rhs.size().first || this->size().second!=rhs.size().second) [[unlikely]] throw std::runtime_error("Wrong vector sizes on copy");
            std::copy(rhs.data.begin(), rhs.data.end(), data.begin());
            return *this;
        }


        
        /*Vector& operator=(const Vector<const dtype>& rhs) requires (!std::is_const_v<dtype>)
        {
            if (this->size().first!=rhs.size().first || this->size().second!=rhs.size().second) [[unlikely]] throw std::runtime_error("Wrong vector sizes on copy");
            std::copy(rhs.data.begin(), rhs.data.end(), data.begin());
            return *this;
        }*/

        Vector(const Vector& rhs):transp(rhs.transp),data(new dtype[rhs.data.size()](0),rhs.data.size()),owns(1)
        {
            *this=rhs;
        }

        Vector(const IMatrix<dtype>& M):owns(1){
            auto p=M.size();
            if (p.second!=1 && p.first!=1) throw std::runtime_error("Too big mstrix to become vector");
            if (p.second==1)
            {
                this->transp=false;
                this->data=std::span<dtype>(new dtype[p.first](0),p.first);
                for (uint32_t i=0;i<p.first;i++)
                    (this->data)[i]=M.gev(i, 0);
            }
            else
            {
                this->transp=true;
                this->data=std::span<dtype>(new dtype[p.second](0),p.second);
                for (uint32_t i=0;i<p.second;i++)
                    (this->data)[i]=M.gev(0, i);
            }
        }

        Vector transposed() const
        {
            Vector res(*this);
            res.transp=1-transp;
            return res;
        }

        dtype lenght() const
        {
            return std::sqrt((*this)|(*this));
        }

        template<typename  T>
        friend T operator*(const Vector<T>&,const Vector<T>&);

        template<typename  T>
        friend T operator|(const Vector<T>&,const Vector<T>&);

        Vector operator+(const Vector &rhs) const
        {
            if (this->transp!=rhs.transp) throw std::runtime_error("Vectors with different state");
            if (this->data.size()!=rhs.data.size()) throw std::runtime_error("Vectors with different size can not be summurized");
            Vector v2(*this);
            for (uint32_t i=0;i<data.size();i++)
                (v2.data)[i]+=(rhs.data)[i];
            return v2;
        }

        Vector operator-(const Vector &rhs) const
        {
            if (this->transp!=rhs.transp) throw std::runtime_error("Vectors with different state");
            if (this->data.size()!=rhs.data.size()) throw std::runtime_error("Vectors with different size can not be summurized");
            Vector v2(*this);
            for (uint32_t i=0;i<data.size();i++)
                (v2.data)[i]-=(rhs.data)[i];
            return v2;
        }

        Vector& operator-=(const Vector &rhs)
        {
            if (this->transp!=rhs.transp) throw std::runtime_error("Vectors with different state");
            if (this->data.size()!=rhs.data.size()) throw std::runtime_error("Vectors with different size can not be summurized");
            for (uint32_t i=0;i<data.size();i++)
                (data)[i]-=(rhs.data)[i];
            return *this;
        }

        Vector operator-() const
        {
            Vector v2(this->size());
            for (uint32_t i=0;i<data.size();i++)
                (v2.data)[i]=-(this->data)[i];
            return v2;
        }

        std::vector<dtype>& get_vector()
        {
            return data;
        }

        dtype& ge(const uint32_t i, const uint32_t j)
        {
            if (!transp){
                //if (j!=0) throw std::runtime_error("Wrong vector usage");
                return (data)[i];
            }
            else {
                //if (i!=0) throw std::runtime_error("Wrong vector usage");
                return (data)[j];
            }
        };

        const dtype& ge(const uint32_t i, const uint32_t j) const
        {
            if (!transp){
                //if (j!=0) throw std::runtime_error("Wrong vector usage");
                return (data)[i];
            }
            else {
                //if (i!=0) throw std::runtime_error("Wrong vector usage");
                return (data)[j];
            }
        };
        
        dtype gev(const uint32_t i, const uint32_t j) const
        {
            if (!transp){
                //if (j!=0) throw std::runtime_error("Wrong vector usage");
                return (data)[i];
            }
            else {
                //if (i!=0) throw std::runtime_error("Wrong vector usage");
                return (data)[j];
            }
        };

        dtype& operator[](uint32_t i)
        {
            return (data)[i];
        }

        const dtype& operator[](uint32_t i) const
        {
            return (data)[i];
        }

        inline const std::pair<uint32_t, uint32_t> size() const {
            if (!transp) return std::make_pair(data.size(), 1);
            else return std::make_pair(1,data.size());
        };

        class Iterator: private IIterator<dtype>
        {
            private:
                uint32_t idx;
            public:
                Iterator(const Vector<dtype>& M, size_t idx_):IIterator<dtype>(M),idx(static_cast<uint32_t>(idx_)){}
                void operator++()
                {
                    idx++;
                    /*if (const Vector<dtype>* sptr = static_cast<const Vector<dtype>*>(this->wptr))
                    {
                        if (idx<sptr->data.size()) 
                    }
                    else
                        throw std::runtime_error("Iterator's object was destroyed");*/
                }

                std::tuple<uint32_t,uint32_t,dtype> operator*()
                {
                    //if (
                    const Vector<dtype>* sptr = static_cast<const Vector<dtype>*>(this->wptr);//)
                    //{
                       // if (idx==sptr->data.size()) throw std::runtime_error("Trying to unname end pointer");
                    return std::make_tuple((sptr->transp?0:idx),
                            (sptr->transp?idx:0),sptr->data[idx]);
                    //}
                    //else
                    //    throw std::runtime_error("Iterator's object was destroyed");
                }

                bool operator!=(const Iterator& it)
                {
                    return (it.idx!=idx)||(this->wptr!=it.wptr);
                }
        };

        Iterator begin()const{return Vector<dtype>::Iterator(*this,0);}
        Iterator end()const{return Vector<dtype>::Iterator(*this,this->data.size());}

        ~Vector()
        {
            if (owns) delete [] data.data();
        }
};

template<typename dtype>
dtype operator*(const Vector<dtype> &lhs,const Vector<dtype> &rhs)
{
    if (!(lhs.transp==1 && rhs.transp==0)) throw std::runtime_error("Vectors multiplication wrong usage");
    if (lhs.data.size()!=rhs.data.size()) throw std::runtime_error("Vectors with different size can not be multiplied");
    dtype res=0;
    for (uint32_t i=0;i<lhs.data.size();i++)
        res+=(lhs.data)[i]*(rhs.data)[i];
    return res;
}

template<typename dtype>
dtype operator|(const Vector<dtype> &lhs,const Vector<dtype> &rhs)
{
    if (!(lhs.transp==rhs.transp)) throw std::runtime_error("Vectors scalar mul wrong usage");
    if (lhs.data.size()!=rhs.data.size()) throw std::runtime_error("Vectors with different size can not be multiplied");
    dtype res=0;
    for (uint32_t i=0;i<lhs.data.size();i++)
        res+=(lhs.data)[i]*(rhs.data)[i];
    return res;
}

template<typename dtype>
dtype norm2(const Vector<dtype> &lhs)
{
    dtype res=0;
    for (uint32_t i=0;i<lhs.size().first || i<lhs.size().second ;i++)
        res=std::max(res,std::abs(lhs[i]));
    return res;
}

};

#endif