#ifndef mdoc_h__
#define mdoc_h__

#include "imatrix.hpp"
#include <cstdint>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

template<typename dtype>
class MDOC: IMatrix<dtype>
{
    private:
        std::map<std::pair<std::uint32_t, std::uint32_t>, dtype> data;
        std::uint32_t n,m;
    public:
        MDOC(){}
        template<IsMatrix<dtype> Mtype>
        MDOC(const Mtype& M):n(M.size().first),m(M.size().second){ //надо переделать с помощью итератора
            for (auto [i,j,c]: M)
                if (c!=0)
                    data[std::make_pair(i, j)]=c;
        }

        MDOC(std::vector<uint32_t> is, std::vector<uint32_t> js,
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
};

#ifndef Dmdoc
extern template class MDOC<double>;
extern template class MDOC<float>;
#endif

#endif