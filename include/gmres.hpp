#ifndef gmres_h__
#define gmres_h__

#include "matrix.hpp"
#include "mvector.hpp"
#include "types.hpp"
#include <cmath>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <utility>


namespace LINARS {

    template<uint32_t M,typename dtype>
    class GMRESdata
    {
        private:
            
        public:
        std::array<std::pair<dtype, dtype>,M> rotations;
        VMatrix<dtype> basis;
        MURtriang<dtype> R;

        GMRESdata(uint32_t v_size):basis(v_size,M),R(M,M){}

        void applyQ(Vector<dtype>& h, const uint32_t& n, const bool& reverse=1)
        {
            dtype a,b;
            if (reverse)
                for (uint32_t i=0;i<n;i++)
                {
                    a=h[i];
                    b=h[i+1];
                    auto [c,s]=this->rotations[i];
                    h[i]=c*a+s*b;
                    h[i+1]=-s*a+c*b;
                }
            else
                for (uint32_t i=n-1;i<n;i--)
                {
                    a=h[i];
                    b=h[i+1];
                    auto [c,s]=this->rotations[i];
                    h[i]=c*a-s*b;
                    h[i+1]=s*a+c*b;
                }
        }
    };

    template<uint32_t M,typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    Vector<dtype> GMRES(const Mtype& A, const Vector<dtype>& b, const Vector<dtype>& x0, uint32_t max_iter=preset_max_iter, dtype max_r=preset_max_r)
    {
        uint32_t n=0;
        GMRESdata<M,dtype> D(x0.size().first);
        Vector<dtype> v=A*x0-b, h(x0.size().first+1), r(x0.size().first+1);
        r[0]=v.lenght();
        D.basis[0]=v/v.lenght();
        for (uint32_t i=1;i<std::min(M,h.size().first);i++)
        {
            //ard
            v=A*D.basis[i-1];
            //std::cout<<i<<" s:\n"<<v;
            for (uint32_t j=0;j<i;j++)
            {
                h[j]=v|D.basis[j];
                v-=h[j]*D.basis[j];
                //std::cout<<i<<" "<<j<<":\n"<<v;
            }
            //std::cout<<"r:\n"<<r<<std::endl;
            h[i]=v.lenght();
            //std::cout<<"vl:"<<h[i]<<" "<<(std::abs(h[i])<max_r)<<" "<<std::abs(h[i])<<" "<<max_r<<std::endl;
            if (std::abs(h[i])>max_r) D.basis[i]=v/h[i];

            //rot and add
            D.applyQ(h, n);
            dtype c=h[i-1]/sqrt(h[i-1]*h[i-1]+h[i]*h[i]),s=-h[i]/sqrt(h[i-1]*h[i-1]+h[i]*h[i]);
            r[i]=std::exchange(r[i-1], r[i-1]*c)*s;
            D.rotations[i-1]=std::make_pair(c, s);
            for (uint32_t j=0;j<i;j++) D.R[j,i-1]=h[j];
            n++;

            //if (std::abs(r[i])<max_r) break;
            //std::cout<<"R:\n"<<D.R<<std::endl<<"r:\n"<<r<<std::endl<<"basis:\n"<<D.basis<<std::endl<<std::endl;
        }
        //std::cout<<D.R<<std::endl<<r<<std::endl<<D.basis<<std::endl;
        for (uint32_t i=n-1;i<n;i--)
        {
            dtype acc=0;
            for (uint32_t j=n-1;j>i;j--)
                acc+=D.R[i,j]*h[j];
            h[i]=(r[i]-acc)/D.R[i,i];
        }
        for (uint32_t i=0;i<v.size().first;i++) v[i]=x0[i];
        //std::cout<<"h:"<<h<<"v1:"<<v<<std::endl;
        for (uint32_t i=0;i<n;i++)
            v-=D.basis[i]*h[i];
        //std::cout<<"v2:"<<v<<std::endl;
        return v;
    }


}


#endif