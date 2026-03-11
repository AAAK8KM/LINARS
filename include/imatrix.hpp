#ifndef imatrix_hpp__
#define imatrix_hpp__

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
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
        
        virtual dtype gev(const uint32_t i, const uint32_t j) const = 0;

        inline virtual const std::pair<uint32_t, uint32_t> size() const = 0;

        virtual Matrix<dtype> operator*(const IMatrix& B) const
        {
            if (this->size().second!=B.size().first) throw std::runtime_error("Matrixes has worng sizes. Can not multiply!");
            Matrix<dtype> C(this->size().first,B.size().second);
            for (uint32_t i=0;i<this->size().first;i++)
                for (uint32_t j=0;j<B.size().second;j++)
                    for (uint32_t k=0;k<this->size().first;k++)
                        C.ge(i, j)+=this->gev(i, k)*B.gev(k, j);
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
Vector<dtype> operator*(const Mtype& A,const Vector<dtype>& B)
{
    if (A.size().second!=B.size().first) throw std::runtime_error("Matrix and vector has worng sizes. Can not multiply!");
    Vector<dtype> C(B.size());
    for (auto [i,j,c]: A)
        C[i]+=c*B[j];
    return C;
}


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

template<typename dtype>
class Vector: public IMatrix<dtype>
{
    private:
        bool transp;
        std::vector<dtype> data;
    public:
        Vector(){}
        Vector(uint32_t n_, bool t=0):transp(t),data(n_,0){}
        Vector(std::pair<uint32_t, uint32_t> s):transp(s.first<s.second),data(s.first<s.second?s.second:s.first)
        {if ((s.first>s.second?s.second:s.first)!=1) throw std::runtime_error("Wrong vector size to create!");}
        Vector(std::vector<dtype>& v, bool t=0):transp(t),data(v){}
        Vector(const IMatrix<dtype>& M){
            auto p=M.size();
            if (p.second!=1 && p.first!=1) throw std::runtime_error("Too big mstrix to become vector");

            if (p.second==1)
            {
                this->transp=false;
                this->data.resize(p.first);
                for (uint32_t i=0;i<p.first;i++)
                    (this->data)[i]=M.gev(i, 0);
            }
            else
            {
                this->transp=true;
                this->data.resize(p.second);
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
class Matrix: public IMatrix<dtype>
{
    private:
        std::vector<dtype> data;
        uint32_t n,m;
    public:
        Matrix(){}
        Matrix(uint32_t n_, uint32_t m_):data(static_cast<std::size_t>(n_)*m_,0),n(n_),m(m_){}
        Matrix(std::pair<uint32_t, uint32_t> s):Matrix(s.first,s.second){}

        template<typename  Mtype>
        requires IsMatrix<dtype, Mtype>
        Matrix(const Mtype& M){
            auto p=M.size();

            this->n=p.first;
            this->m=p.second;

            this->data.resize(m*n);

            for (auto [i,j,c] : M)
                    this->data[static_cast<std::size_t>(i)*p.second+j]=c;
                
        }

        dtype& ge(const uint32_t i, const uint32_t j)
        {
            //std::cout<<static_cast<std::size_t>(i)<<" "<<j<<std::endl;
            return data[static_cast<std::size_t>(i)*m+j];
        };

        const dtype& ge(const uint32_t i, const uint32_t j) const
        {
            return data[static_cast<std::size_t>(i)*m+j];
        };
        
        dtype gev(const uint32_t i, const uint32_t j) const
        {
            return data[static_cast<std::size_t>(i)*m+j];
        };

        inline const std::pair<uint32_t, uint32_t> size() const {return std::make_pair(n, m);};


        class Iterator: private IIterator<dtype>
        {
            private:
                size_t idx;
                uint32_t it,jt;
            public:
                Iterator(const Matrix<dtype>& M, size_t idx_):IIterator<dtype>(M),idx(idx_),
                        it(static_cast<uint32_t>(idx_/M.m)),jt(static_cast<uint32_t>(idx_/M.m)){}
                void operator++()
                {
                    const Matrix<dtype>* sptr = static_cast<const Matrix<dtype>*>(this->wptr);
                    idx++;
                    jt++;
                    if (jt==sptr->m)
                    {
                        jt=0;
                        it++;
                    }
                    /*if (const Matrix<dtype>* sptr = static_cast<const Matrix<dtype>*>(this->wptr))
                    {
                        if (idx<sptr->data.size()) 
                    }
                    else
                        throw std::runtime_error("Iterator's object was destroyed");*/
                }

                std::tuple<uint32_t,uint32_t,dtype> operator*()
                {
                    const Matrix<dtype>* sptr = static_cast<const Matrix<dtype>*>(this->wptr);
                        //if (idx==sptr->data.size()) throw std::runtime_error("Trying to unname end pointer");
                    //std::size_t p=idx/sptr->m;
                    return std::make_tuple(static_cast<uint32_t>(it),
                            static_cast<uint32_t>(jt),sptr->data[idx]);
                }

                bool operator!=(const Iterator& itr)
                {
                    return (itr.idx!=idx)||(this->wptr!=itr.wptr);
                }
        };

        Iterator begin()const{return Matrix<dtype>::Iterator(*this,0);}
        Iterator end()const{return Matrix<dtype>::Iterator(*this,this->data.size());}
};

template<typename dtype>
class VMatrix: public IMatrix<dtype>
{
    private:
        std::vector<Vector<dtype>> data;
        uint32_t n,m;
    public:
        VMatrix(){}
        VMatrix(uint32_t n_, uint32_t m_):data(static_cast<std::size_t>(m_),Vector<dtype>(n_)),n(n_),m(m_){}
        VMatrix(std::pair<uint32_t, uint32_t> s):VMatrix(s.first,s.second){}

        template<typename  Mtype>
        requires IsMatrix<dtype, Mtype>
        VMatrix(const Mtype& M){
            auto p=M.size();

            this->n=p.first;
            this->m=p.second;

            this->data.resize(m);
            for (uint32_t i=0;i<m;i++)
                data[i]=Vector<dtype>(n);

            for (auto [i,j,c] : M)
                    this->data[j][i]=c;
                
        }

        dtype& ge(const uint32_t i, const uint32_t j)
        {
            //std::cout<<static_cast<std::size_t>(i)<<" "<<j<<std::endl;
            return data[j][i];
        };

        const dtype& ge(const uint32_t i, const uint32_t j) const
        {
            return data[j][i];
        };
        
        dtype gev(const uint32_t i, const uint32_t j) const
        {
            return data[j][i];
        };

        Vector<dtype>& gv(const uint32_t j)
        {
            return data[j];
        };

        Vector<dtype>& operator[](const uint32_t j)
        {
            return data[j];
        };

        const Vector<dtype>& operator[](const uint32_t j) const
        {
            return data[j];
        };

        inline const std::pair<uint32_t, uint32_t> size() const {return std::make_pair(n, m);};


        class Iterator: private IIterator<dtype>
        {
            private:
                size_t idx;
                uint32_t it,jt;
            public:
                Iterator(const VMatrix<dtype>& M, size_t idx_):IIterator<dtype>(M),idx(idx_),
                        it(static_cast<uint32_t>(idx_/M.m)),jt(static_cast<uint32_t>(idx_/M.m)){}
                void operator++()
                {
                    const VMatrix<dtype>* sptr = static_cast<const VMatrix<dtype>*>(this->wptr);
                    idx++;
                    it++;
                    if (it==sptr->n)
                    {
                        it=0;
                        jt++;
                    }
                    /*if (const Matrix<dtype>* sptr = static_cast<const Matrix<dtype>*>(this->wptr))
                    {
                        if (idx<sptr->data.size()) 
                    }
                    else
                        throw std::runtime_error("Iterator's object was destroyed");*/
                }

                std::tuple<uint32_t,uint32_t,dtype> operator*()
                {
                    const VMatrix<dtype>* sptr = static_cast<const VMatrix<dtype>*>(this->wptr);
                        //if (idx==sptr->data.size()) throw std::runtime_error("Trying to unname end pointer");
                    //std::size_t p=idx/sptr->m;
                    return std::make_tuple(static_cast<uint32_t>(it),
                            static_cast<uint32_t>(jt),sptr->data[jt][it]);
                }

                bool operator!=(const Iterator& itr)
                {
                    return (itr.idx!=idx)||(this->wptr!=itr.wptr);
                }
        };

        Iterator begin()const{return VMatrix<dtype>::Iterator(*this,0);}
        Iterator end()const{return VMatrix<dtype>::Iterator(*this,this->data.size()*n);}
};

#ifndef Dmatrix
extern template class Matrix<double>;
extern template class Matrix<float>;
extern template class VMatrix<double>;
extern template class VMatrix<float>;
#endif

#ifndef Dvector
extern template class Vector<double>;
extern template class Vector<float>;
#endif

}

#endif