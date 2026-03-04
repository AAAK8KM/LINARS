#ifndef qrdec_h__
#define qrdec_h__

#include "imatrix.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <vector>
#include <ranges>

namespace LINARS {

template<typename dtype>
class HausholdOP
{
    private:
        uint32_t start;
        Vector<dtype> data;
    public:
        HausholdOP(Vector<dtype> v, uint32_t start_):start(start_),data(v){};

        uint32_t get_start() {return start;};
        Vector<dtype> operator()(const Vector<dtype>& v)
        {
            return v-(2*(data.transposed_shared()*v)/data.lenght())*data;
        }

        VMatrix<dtype> operator()(const VMatrix<dtype>& M)
        {
            VMatrix<dtype> res(M);
            for (uint32_t i=start;i<res.size().second;i++) res[i]=(*this)(res[i]);
            return res;
        }
};

template<typename dtype>
class QHaushold
{
    private:
    bool order=1;
    std::vector<HausholdOP<dtype>> data;
    public:
    VMatrix<dtype> operator()(const VMatrix<dtype>& M)
    {
        VMatrix<dtype> res(M);
        for (uint32_t idx=(order?0:data.size()-1);idx<data.size();(order?idx++:idx--))
            for (uint32_t i=data[idx].get_start();i<M.size().second;i++) M[i]=data[idx](M[i]);
        return res;
    }
    bool get_order(){return order;}
    void reverse(){order=!order;}

    template<typename  T>
    friend std::pair<QHaushold<T>, VMatrix<T>> QRdecomposition(const IMatrix<T>& M);
};


}

#endif