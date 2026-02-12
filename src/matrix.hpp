#ifndef matrix_hpp__
#define matrix_hpp__


#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

template<typename dtype>
class Matrix;


template<typename dtype>
class Imatrix
{
    public:
        virtual dtype& ge(const std::size_t i, const std::size_t j) = 0;
        virtual const dtype& ge(const std::size_t i, const std::size_t j) const = 0;
        
        virtual dtype gev(const std::size_t i, const std::size_t j) const = 0;

        inline virtual const std::pair<std::size_t, std::size_t> size() const = 0;
        //virtual dtype det();//later

        Matrix<dtype> operator*(const Imatrix& B) const
        {
            if (this->size().second!=B.size().first) throw std::runtime_error("Matrixes has worng sizes. Can not multiply!");
            Matrix<dtype> C(this->size().first,B.size().second);
            for (std::size_t i=0;i<this->size().first;i++)
                for (std::size_t j=0;j<B.size().second;j++)
                    for (std::size_t k=0;k<this->size().first;k++)
                        C.ge(i, j)+=this->gev(i, k)*B.gev(k, j);
            return C;
        }

        bool operator==(const Imatrix& B) const
        {
            if (this->size().second!=B.size().second) return false;
            if (this->size().first!=B.size().first) return false;

            for (std::size_t i=0;i<this->size().first;i++)
                for (std::size_t j=0;j<this->size().second;j++)
                    if (std::abs(this->gev(i,j)-B.gev(i, j))>1e-9)
                        return false;
            return true;
        }

        bool operator!=(const Imatrix& B) const
        {
            return !(*this==B);
        }
};

// To do later
/*template<typename dtype, std::size_t N, std::size_t M>
class Imatrixf
{
    inline const std::pair<std::size_t, std::size_t> size() {return std::make_pair(N, M);}
};*/

template<typename dtype>
class Matrix: public Imatrix<dtype>
{
    private:
        std::vector<dtype> data;
        std::uint32_t n,m;
    public:
        Matrix(){}
        Matrix(std::uint32_t n, std::uint32_t m):data(n*m,0),n(n),m(m){}
        Matrix(std::pair<std::uint32_t, std::uint32_t> s):Matrix(s.first,s.second){}
        Matrix(Imatrix<dtype>& M){
            auto p=M.size();
            
            if (std::max(p.first,p.second)>UINT32_MAX)  throw std::runtime_error("Too big matrix!");
            
            this->n=p.first;
            this->m=p.second;

            this->data.resize(m*n);

            for (std::size_t i=0;i<p.first;i++)
                for (std::size_t j=0;j<p.second;j++)
                    this->data[i*p.second+j]=M.gev(i, j);
                
        }

        Matrix(Matrix& M):data(M.data),n(M.n),m(M.m){}
        Matrix& operator=(const Matrix& M)
        {
            this->data=M.data;
            this->n=M.n;
            this->m=M.m;
            return *this;
        }
        Matrix(Matrix&& M)
        {
            this->data=std::move(M.data);
            this->n=std::move(M.n);
            this->m=std::move(M.m);
        };
        Matrix& operator=(const Matrix&& M)
        {
            this->data=std::move(M.data);
            this->n=std::move(M.n);
            this->m=std::move(M.m);
            return *this;
        }

        dtype& ge(const std::size_t i, const std::size_t j)
        {
            return data[i*m+j];
        };

        const dtype& ge(const std::size_t i, const std::size_t j) const
                {
            return data[i*m+j];
        };
        
        dtype gev(const std::size_t i, const std::size_t j) const
        {
            return data[i*m+j];
        };

        inline const std::pair<std::size_t, std::size_t> size() const {return std::make_pair(n, m);};

        ~Matrix(){};
};


template<typename dtype>
class M3diag: public Imatrix<dtype>
{
    private:
        std::vector<dtype> data;
    public:
        M3diag(){}
        M3diag(std::size_t n):data(3*n,0){}
        M3diag(Imatrix<dtype>& m):data(3*m.size().first,0){
            auto p=m.size();
            if (p.first!=p.second) throw std::runtime_error("Can not create 3 diag matrix from not square!");
            for (std::size_t i=0;i<p.first;i++)
                for (std::size_t j=0;j<p.first;j++)
                {
                    if ((i > j ? i - j : j - i) > 1)
                    {
                        if (m.gev(i, j)){
                            throw std::runtime_error("Not valid matrix!");
                        }
                    }
                    else
                    {
                        if (i-j==1)
                            this->data[i*3-1]=m.gev(i, j);
                        if (i==j)
                            this->data[i*3]=m.gev(i, j);
                        if (j-i==1)
                            this->data[i*3+1]=m.gev(i, j);
                    }
                }
        }
        M3diag(M3diag& m):data(m.data){}
        M3diag& operator=(const M3diag& m)
        {
            this->data=m.data;
            return *this;
        }
        M3diag(M3diag&& m)
        {
            this->data=std::move(m.data);
        };
        M3diag& operator=(const M3diag&& m)
        {
            this->data=std::move(m.data);
            return *this;
        }

        dtype& ge(const std::size_t i, const std::size_t j)
        {
            if ((i > j ? i - j : j - i) > 1)
            {
                std::cout<<i<<" "<<j<<std::endl;
                throw std::runtime_error("Not valid element!");
            }
            return data[2*i+j];
            //else return data[3*i+1];
        };

        const dtype& ge(const std::size_t i, const std::size_t j) const
                {
            if ((i > j ? i - j : j - i) > 1)
            {
                std::cout<<i<<" "<<j<<std::endl;
                throw std::runtime_error("Not valid element!.");
            }
            return data[2*i+j];
        };
        
        dtype gev(const std::size_t i, const std::size_t j) const
        {
            if ((i > j ? i - j : j - i) > 1) return 0;
            //if (i>j) return data[4*i-j];
            return data[2*i+j];
        };

        inline const std::pair<std::size_t, std::size_t> size() const {return std::make_pair(data.size()/3, data.size()/3);};

        void flush()
        {
            for (auto d :data)
            std::cout<<d<<" ";
            std::cout<<std::endl;
        }
        ~M3diag(){};
};


// To do later
/*template<typename dtype, std::size_t N>
class M3diagf: public Imatrixf<dtype, N, N>
{
    private:
        dtype data[3*N];
};*/

#ifndef Dmatrix
extern template class Matrix<double>;
extern template class Matrix<float>;

extern template class M3diag<double>;
extern template class M3diag<float>;
#endif

#endif