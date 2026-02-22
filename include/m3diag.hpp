#ifndef m3diag_hpp__
#define m3diag_hpp__

#include "imatrix.hpp"

template<typename dtype>
class M3diag: public Imatrix<dtype>
{
    private:
        std::vector<dtype> data;
    public:
        M3diag(){}
        M3diag(std::uint32_t n):data(3*n,0){}
        M3diag(Imatrix<dtype>& m):data(3*m.size().first,0){
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
        M3diag(M3diag& m):data(m.data){}

        dtype& ge(const std::uint32_t i, const std::uint32_t j)
        {
            if ((i > j ? i - j : j - i) > 1)
            {
                std::cout<<i<<" "<<j<<std::endl;
                throw std::runtime_error("Not valid element!");
            }
            return data[2*static_cast<std::size_t>(i)+j];
            //else return data[3*i+1];
        };

        const dtype& ge(const std::uint32_t i, const std::uint32_t j) const
                {
            if ((i > j ? i - j : j - i) > 1)
            {
                std::cout<<i<<" "<<j<<std::endl;
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

        inline const std::pair<std::uint32_t, std::uint32_t> size() const {return std::make_pair(data.size()/3, data.size()/3);};

        ~M3diag(){};
};

#ifndef Dmatrix
extern template class M3diag<double>;
extern template class M3diag<float>;
#endif


#endif