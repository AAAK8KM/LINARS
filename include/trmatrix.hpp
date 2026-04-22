#ifndef trmatrix_hpp__
#define trmatrix_hpp__


#include "imatrix.hpp"
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace LINARS {

template <typename dtype>
class MURtriang: public IMatrix<dtype>
{
    private:
        uint32_t n,m;
        std::vector<dtype> data;

    public:
        MURtriang(const uint32_t& n, const uint32_t& m):MURtriang(std::make_pair(n, m)){}

        MURtriang(const std::pair<uint32_t, uint32_t>& size):
        n(size.first),m(size.second),data(std::min(n,m)*(2*m-std::min(n,m)+1)/2){}

        template<typename  Mtype>
        requires IsMatrix<dtype, Mtype>
        MURtriang(const Mtype& M):MURtriang(M.size())
        {
            for (auto [i,j,c] : M)
            {
                if (c!=0)
                {
                    if (i>j) throw std::invalid_argument("Wrong matrix to make threediag");
                    else this->ge(i, j)=c;
                }
            }
        }

        dtype& ge(const std::uint32_t i, const std::uint32_t j)
        {
            return data[(j*(j+1))/2+i];
        };

        const dtype& ge(const std::uint32_t i, const std::uint32_t j) const
        {
            return data[(j*(j+1))/2+i];
        };
        
        dtype gev(const std::uint32_t i, const std::uint32_t j) const
        {
            if (i<=j && j<m) return ge(i, j);
            else return 0;
        };

        const std::pair<uint32_t, uint32_t> size() const {return  std::make_pair(n, m);}

        class Iterator: private IIterator<dtype>
        {
            private:
                uint32_t idx;
                uint32_t i,j;
            public:
                Iterator(const MURtriang<dtype>& M, uint32_t idx_):IIterator<dtype>(M),idx(idx_),i(0),j(0){}
                void operator++()
                {
                    //if (
                    const MURtriang<dtype>* sptr = static_cast<const MURtriang<dtype>*>(this->wptr);//)
                    //{
                        /*if (idx!=sptr->vals.size()) 
                        {*/
                            idx++;
                            if (i==j)
                            {
                                j++;
                                i=0;
                            }
                            else i++;
                        //}
                    //}
                    //else
                    //    throw std::runtime_error("Iterator's object was destroyed");
                }

                std::tuple<uint32_t,uint32_t,dtype> operator*()
                {
                    //if (
                    const MURtriang<dtype>* sptr = static_cast<const MURtriang<dtype>*>(this->wptr);//)
                    //{
                       // if (idx==sptr->vals.size()) throw std::runtime_error("Trying to unname end pointer");
                        return std::make_tuple(i,j,sptr->data[idx]);
                    //}
                    //else
                    //    throw std::runtime_error("Iterator's object was destroyed");
                }

                bool operator!=(const Iterator& rhs)
                {
                    return (idx!=rhs.idx)||(this->wptr!=rhs.wptr);
                }
        };

        Iterator begin()const{return MURtriang<dtype>::Iterator(*this,0);}
        Iterator end()const{return MURtriang<dtype>::Iterator(*this,static_cast<uint32_t>(data.size()));}
};


}

#endif