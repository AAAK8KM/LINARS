#ifndef mcsr_h__
#define mcsr_h__


#include "imatrix.hpp"
#include <algorithm>
#include <cstdint>
#include <tuple>
#include <utility>
#include <vector>

namespace LINARS {

template <typename dtype>
class MCSR: public IMatrix<dtype>
{
    private:
        std::vector<std::pair<uint32_t,dtype>> vals;
        std::vector<uint32_t> rows;
        uint32_t n,m;

        static bool comp(const std::pair<uint32_t, dtype>& a,const std::pair<uint32_t, dtype>& b)
        {
            return  (a.first<b.first);
        }
    public:
        MCSR(const std::pair<uint32_t, uint32_t>& size):n(size.first),m(size.second){}

        template<typename  Mtype>
        requires IsMatrix<dtype, Mtype>
        MCSR(const Mtype& M):rows(M.size().first+1,0),n(rows.size()-1),m(M.size().second)
        {
            std::vector<std::tuple<uint32_t,uint32_t,dtype>> v;
            for (auto [i,j,c] : M)
            {
                if (c!=0)
                {
                    rows[i+1]+=1;
                    v.push_back(std::make_tuple(i,j,c));
                }
            }
            std::sort(v.begin(),v.end());
            vals.resize(v.size());
            for (uint32_t i=0;i<v.size();i++)
                vals[i]=std::make_pair(std::get<1>(v[i]), std::get<2>(v[i]));
            /*for (uint32_t i=0;i<=n;i++)
                std::cout<<rows[i]<<" ";
            std::cout<<std::endl;*/
            for (uint32_t i=1;i<=n;i++)
                rows[i]+=rows[i-1];
            /*for (uint32_t i=0;i<=n;i++)
                std::cout<<rows[i]<<" ";
            std::cout<<std::endl;
            for (uint32_t i=0;i<vals.size();i++)
            {
                std::cout<<"("<<vals[i].first<<" "<<vals[i].second<<") ";
            }
            std::cout<<std::endl;*/
        }

        dtype& ge(const std::uint32_t i, const std::uint32_t j)
        {
            auto it=std::lower_bound(vals.begin()+rows[i],vals.begin()+rows[i+1],std::make_pair(j, 0),comp);
            if (j==it->first && it!=vals.begin()+rows[i+1]) return it->second;
            throw std::runtime_error("Not valid element!");
        };

        const dtype& ge(const std::uint32_t i, const std::uint32_t j) const
        {
            auto it=std::lower_bound(vals.begin()+rows[i],vals.begin()+rows[i+1],std::make_pair(j, 0),comp);
            if (j==it->first && it!=vals.begin()+rows[i+1]) return it->second;
            throw std::runtime_error("Not valid element!");
        };
        
        dtype gev(const std::uint32_t i, const std::uint32_t j) const
        {
            auto it=std::lower_bound(vals.begin()+rows[i],vals.begin()+rows[i+1],std::make_pair(j, 0),comp);
            if (j==it->first && it!=vals.begin()+rows[i+1]) return it->second;
            else return 0;
        };

        const std::pair<uint32_t, uint32_t> size() const {return  std::make_pair(n, m);}

        class Iterator: private IIterator<dtype>
        {
            private:
                uint32_t idx;
                uint32_t i;
            public:
                Iterator(const MCSR<dtype>& M, uint32_t idx_):IIterator<dtype>(M),idx(idx_),i(0){while (M.rows[i+1]<=idx && (i+1)!=M.rows.size()-1) i++;}
                void operator++()
                {
                    //if (
                    const MCSR<dtype>* sptr = static_cast<const MCSR<dtype>*>(this->wptr);//)
                    //{
                        /*if (idx!=sptr->vals.size()) 
                        {*/
                            idx++;
                            while (sptr->rows[i+1]<=idx && (i+1)!=sptr->rows.size()-1) i++;
                        //}
                    //}
                    //else
                    //    throw std::runtime_error("Iterator's object was destroyed");
                }

                std::tuple<uint32_t,uint32_t,dtype> operator*()
                {
                    //if (
                    const MCSR<dtype>* sptr = static_cast<const MCSR<dtype>*>(this->wptr);//)
                    //{
                       // if (idx==sptr->vals.size()) throw std::runtime_error("Trying to unname end pointer");
                        auto p=sptr->vals[idx];
                        return std::make_tuple(i,p.first,p.second);
                    //}
                    //else
                    //    throw std::runtime_error("Iterator's object was destroyed");
                }

                bool operator!=(const Iterator& rhs)
                {
                    return (idx!=rhs.idx)||(this->wptr!=rhs.wptr);
                }
        };

        Iterator begin()const{return MCSR<dtype>::Iterator(*this,0);}
        Iterator end()const{return MCSR<dtype>::Iterator(*this,static_cast<uint32_t>(vals.size()));}
};

#ifndef Dmcsr
extern template class MCSR<double>;
extern template class MCSR<float>;
#endif

}

#endif