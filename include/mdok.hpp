#ifndef mdok_h__
#define mdok_h__

#include "imatrix.hpp"
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace LINARS {

template<typename dtype>
class MDOK: public IMatrix<dtype>
{
    private:
        std::map<std::pair<std::uint32_t, std::uint32_t>, dtype> data;
        std::uint32_t n,m;
    public:
        MDOK(const std::pair<uint32_t, uint32_t>& size):n(size.first),m(size.second){}
        
        template<typename  Mtype>
        requires IsMatrix<dtype, Mtype>
        MDOK(const Mtype& M):n(M.size().first),m(M.size().second){ //надо переделать с помощью итератора
            for (auto [i,j,c]: M)
                if (c!=0)
                    data[std::make_pair(i, j)]=c;
        }

        MDOK(std::vector<uint32_t> is, std::vector<uint32_t> js,
             std::vector<dtype> cs, uint32_t n_, uint32_t m_):n(n_),m(m_)
        {
            if (is.size()!=js.size() || is.size()!=cs.size()) throw std::runtime_error("bad input data");
            for (uint32_t i=0;i<is.size();i++)
                data[std::make_pair(is[i], js[i])]=cs[i];
        }

        void set(uint32_t i, uint32_t j, dtype c)
        {
            if (i>=n || j>=m) throw std::out_of_range("Invalid indexes");
            data[std::make_pair(i, j)]=c;
        }
        
        dtype& ge(const std::uint32_t i, const std::uint32_t j)
        {
            return data[std::make_pair(i, j)];
        };

        const dtype& ge(const std::uint32_t i, const std::uint32_t j) const
        {
            return data.at(std::make_pair(i, j));
        };
        
        dtype gev(const std::uint32_t i, const std::uint32_t j) const
        {
            auto p= data.find(std::make_pair(i, j));
            if (p!=data.end())
                return p->second;
            else
                return 0;
        };

        const std::pair<uint32_t, uint32_t> size() const {return  std::make_pair(n, m);}

        class Iterator: private IIterator<dtype>
        {
            private:
                std::map<std::pair<std::uint32_t, std::uint32_t>, dtype>::const_iterator it;
            public:
                Iterator(const MDOK<dtype>& M, bool e):IIterator<dtype>(M),it(!e?M.data.begin():M.data.end()){}
                void operator++()
                {
                    //if (
                    //const MDOK<dtype>* sptr = static_cast<const MDOK<dtype>*>(this->wptr))
                    //{
                    //    if (it!=sptr->data.end()) 
                    it++;
                    //}
                    //else
                    //    throw std::runtime_error("Iterator's object was destroyed");
                }

                std::tuple<uint32_t,uint32_t,dtype> operator*()
                {
                    //if (const MDOK<dtype>* sptr = static_cast<const MDOK<dtype>*>(this->wptr))
                    //{
                    //    if (it==sptr->data.end()) throw std::runtime_error("Trying to unname end pointer");
                        return std::make_tuple(it->first.first,it->first.second,it->second);
                    //}
                    //else
                    //    throw std::runtime_error("Iterator's object was destroyed");
                }

                bool operator!=(const Iterator& rhs)
                {
                    return (it!=rhs.it)||(this->wptr!=rhs.wptr);
                }
        };

        Iterator begin()const{return MDOK<dtype>::Iterator(*this,0);}
        Iterator end()const{return MDOK<dtype>::Iterator(*this,1);}
};

#ifndef Dmdok
extern template class MDOK<double>;
extern template class MDOK<float>;
#endif

}

#endif