#ifndef matrix_hpp__
#define matrix_hpp__

#include "imatrix.hpp"
#include "mvector.hpp"
#include <vector>

namespace LINARS {


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
                        it(static_cast<uint32_t>(idx_/M.m)),jt(static_cast<uint32_t>(idx_%M.m)){}
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
        std::vector<dtype> data;
        uint32_t n,m;
    public:
        VMatrix(){}
        VMatrix(uint32_t n_, uint32_t m_):data(n_*m_),n(n_),m(m_){}
        VMatrix(std::pair<uint32_t, uint32_t> s):VMatrix(s.first,s.second){}

        VMatrix& operator=(const VMatrix& rhs) {
            if (this->size().first!=rhs.size().first || this->size().second!=rhs.size().second) [[unlikely]] throw std::runtime_error("Wrong VMatrix sizes on copy");
            std::copy(rhs.data.begin(), rhs.data.end(), data.begin());
            return *this;
        }

        VMatrix(const VMatrix& rhs) = default;

        template<typename  Mtype>
        requires IsMatrix<dtype, Mtype>
        VMatrix(const Mtype& M){
            auto p=M.size();

            this->n=p.first;
            this->m=p.second;

            this->data.resize(m*n);

            for (auto [i,j,c] : M)
                    this->data[j*n+i]=c;
                
        }

        dtype& ge(const uint32_t i, const uint32_t j)
        {
            //std::cout<<static_cast<std::size_t>(i)<<" "<<j<<std::endl;
            return data[j*n+i];
        };

        const dtype& ge(const uint32_t i, const uint32_t j) const
        {
            return data[j*n+i];
        };
        
        dtype gev(const uint32_t i, const uint32_t j) const
        {
            return data[j*n+i];
        };

        using IMatrix<dtype>::operator[];

        VMatrix operator-(const VMatrix<dtype> rhs) const
        {
            VMatrix<dtype> res(*this);
            for (uint32_t i=0;i<this->m;i++)
                res[i]=(*this)[i]-rhs[i];
            return res;
        }

        Vector<dtype> operator[](const uint32_t j)
        {
            return Vector<dtype>(data.data()+j*n,n);
        };

        const Vector<dtype> operator[](const uint32_t j) const
        {
            return Vector<dtype>(data.data()+j*n,n);
        };

        inline const std::pair<uint32_t, uint32_t> size() const {return std::make_pair(n, m);};

        class Iterator: private IIterator<dtype>
        {
            private:
                size_t idx;
                uint32_t it,jt;
            public:
                Iterator(const VMatrix<dtype>& M, size_t idx_):IIterator<dtype>(M),idx(idx_),
                        it(static_cast<uint32_t>(idx_%M.n)),jt(static_cast<uint32_t>(idx_/M.n)){}
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
                            static_cast<uint32_t>(jt),sptr->data[jt*sptr->n+it]);
                }

                bool operator!=(const Iterator& itr)
                {
                    return (itr.idx!=idx)||(this->wptr!=itr.wptr);
                }
        };

        Iterator begin()const{return VMatrix<dtype>::Iterator(*this,0);}
        Iterator end()const{return VMatrix<dtype>::Iterator(*this,this->data.size());}
};

}

#endif