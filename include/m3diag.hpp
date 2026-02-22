#ifndef m3diag_hpp__
#define m3diag_hpp__

#include "imatrix.hpp"
#include <cstddef>
#include <cstdint>

template<typename dtype>
class M3diag: public IMatrix<dtype>
{
    private:
        std::vector<dtype> data;
    public:
        M3diag(){}
        M3diag(std::uint32_t n):data(3*n-2,0){}
        M3diag(IMatrix<dtype>& m):data(3*m.size().first-2,0){
            auto p=m.size();
            if (p.first!=p.second) throw std::runtime_error("Can not create 3 diag matrix from not square!");
            for (std::uint32_t i=0;i<p.first;i++)
                for (std::uint32_t j=0;j<p.first;j++)
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
                            this->data[static_cast<std::size_t>(i)*3-1]=m.gev(i, j);
                        if (i==j)
                            this->data[static_cast<std::size_t>(i)*3]=m.gev(i, j);
                        if (j-i==1)
                            this->data[static_cast<std::size_t>(i)*3+1]=m.gev(i, j);
                    }
                }
        }

        dtype& ge(const std::uint32_t i, const std::uint32_t j)
        {
            if ((i > j ? i - j : j - i) > 1)
            {
                //std::cout<<i<<" "<<j<<std::endl;
                throw std::runtime_error("Not valid element!");
            }
            return data[2*static_cast<std::size_t>(i)+j];
            //else return data[3*i+1];
        };

        const dtype& ge(const std::uint32_t i, const std::uint32_t j) const
                {
            if ((i > j ? i - j : j - i) > 1)
            {
                //std::cout<<i<<" "<<j<<std::endl;
                throw std::runtime_error("Not valid element!.");
            }
            return data[2*static_cast<std::size_t>(i)+j];
        };
        
        dtype gev(const std::uint32_t i, const std::uint32_t j) const
        {
            if ((i > j ? i - j : j - i) > 1) return 0;
            //if (i>j) return data[4*i-j];
            return data[2*static_cast<std::size_t>(i)+j];
        };

        inline const std::pair<std::uint32_t, std::uint32_t> size() const {return std::make_pair((data.size()+2)/3, (data.size()+2)/3);};


        class Iterator: private IIterator<dtype>
        {
            private:
                uint32_t idx;
            public:
                Iterator(const M3diag<dtype>& M, size_t idx_):IIterator<dtype>(M),idx(static_cast<uint32_t>(idx_)){}
                void operator++()
                {
                    if (const M3diag<dtype>* sptr = dynamic_cast<const M3diag<dtype>*>(this->wptr))
                    {
                        if (idx<sptr->data.size()) idx++;
                    }
                    else
                        throw std::runtime_error("Iterator's object was destroyed");
                }

                std::tuple<uint32_t,uint32_t,dtype> operator*()
                {
                    if (const M3diag<dtype>* sptr = dynamic_cast<const M3diag<dtype>*>(this->wptr))
                    {
                        if (idx==sptr->data.size()) throw std::runtime_error("Trying to unname end pointer");
                        return std::make_tuple((idx+1)/3,
                            (idx+1)/3+(idx+1)%3-1,sptr->data[idx]);
                    }
                    else
                        throw std::runtime_error("Iterator's object was destroyed");
                }

                bool operator!=(const Iterator& it)
                {
                    return (it.idx!=idx)||(this->wptr!=it.wptr);
                }
        };

        Iterator begin(){return M3diag<dtype>::Iterator(*this,0);}
        Iterator end(){return M3diag<dtype>::Iterator(*this,this->data.size());}
        
};

#ifndef Dm3diag
extern template class M3diag<double>;
extern template class M3diag<float>;
#endif


#endif