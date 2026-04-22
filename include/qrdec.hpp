#ifndef qrdec_h__
#define qrdec_h__

#include "matrixes.hpp"
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>
//#include "t2m.hpp"

namespace LINARS {

template<typename dtype>
class HausholdOP
{
    private:
        uint32_t start;
        const Vector<dtype> data;
        const dtype l;
    public:
        HausholdOP(const Vector<dtype>& v, uint32_t start_):start(start_),data(v),l((v|v)/2){};

        uint32_t get_start() {return start;};
        Vector<dtype> operator()(const Vector<dtype>& v)
        {
            //std::cout<<data<<std::endl;
            if (v.size().first>v.size().second)
                return v-((data|v)/l)*data;
            else
                return v-((v*data)/l)*data.transposed();
        }

        VMatrix<dtype> operator()(const VMatrix<dtype>& M)
        {
            VMatrix<dtype> res(M);
            for (uint32_t i=0;i<res.size().second;i++) res[i]=(*this)(res[i]);
            return res;
        }
};

template<typename dtype>
class QHaushold
{
    private:
    bool order=0;
    std::vector<HausholdOP<dtype>> data;
    public:


    VMatrix<dtype> opt(const VMatrix<dtype>& M)
    {
        VMatrix<dtype> res(M);
        for (uint32_t idx=(order?0:data.size()-1);idx<data.size();(order?idx++:idx--))
            for (uint32_t i=data[idx].get_start();i<M.size().second;i++) res[i]=data[idx](res[i]);
        return res;
    }

    template <typename  Mtype>
    requires IsMatrix<dtype, Mtype> && (!std::is_same<Mtype, Vector<dtype>>::value)
    VMatrix<dtype> operator()(const Mtype& M)
    {
        VMatrix<dtype> res(M);
        for (uint32_t idx=(order?0:data.size()-1);idx<data.size();(order?idx++:idx--))
            for (uint32_t i=0;i<M.size().second;i++) res[i]=data[idx](res[i]);
        return res;
    }

    Vector<dtype> operator()(const Vector<dtype>& v)
    {
        Vector<dtype> res(v);
        for (uint32_t idx=(order?0:data.size()-1);idx<data.size();(order?idx++:idx--))
            res=data[idx](res);
        return res;
    }

    bool get_order(){return order;}
    void reverse(){order=!order;}

    template <typename T,typename  Mtype>
    requires IsMatrix<T, Mtype>
    friend std::pair<QHaushold<T>, VMatrix<T>> QRdecompositionH(const Mtype& M);
};

template <typename dtype,typename  Mtype>
requires IsMatrix<dtype, Mtype>
std::pair<QHaushold<dtype>, VMatrix<dtype>> QRdecompositionH(const Mtype& M)
{
    QHaushold<dtype> Q;
    VMatrix<dtype> R(M);
    Vector<dtype> v(M.size().first);

    if (R.size().first<R.size().second) throw std::runtime_error("For QR decomposition matrix with sizes MxN is needed to be M>N");
    for (uint32_t i=0;i<R.size().second-1;i++)
    {
        //std::cout<<R<<std::endl;
        v=R[i];
        for (uint32_t idx=0;idx<i;idx++)
            v[idx]=0;
        v[i]+=v.lenght()*(v[i]>0?1:-1);
        if (v.lenght()==0) continue;
        //std::cout<<v<<std::endl;
        Q.data.emplace_back(v,i);
        for (uint32_t idx=i;idx<R.size().second;idx++)
            R[idx]=Q.data.back()(R[idx]);
    }
    return std::make_pair(Q, R);
}

}

#endif