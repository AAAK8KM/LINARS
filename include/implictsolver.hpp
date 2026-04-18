#ifndef implictsolver_h__
#define implictsolver_h__

#include "matrixes.hpp"
#include "t2m.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <chebishev.hpp>
#include <ostream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace LINARS {

    constexpr const uint32_t preset_max_iter=10000;
    constexpr const long double preset_max_r=1e-9;

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> JakobiSolver(const Mtype& A, const VMatrix<dtype>& b, uint32_t max_iter=preset_max_iter, dtype max_r=preset_max_r)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();
        uint32_t iter=0;
        Vector<dtype> x(b.size().first);
        Vector<dtype> r(b.size().first);
        while (iter++<max_iter && mes_r>max_r) {
            mes_r=0;
            for (uint32_t si=0;si<b.size().second;si++)
            {
                x=b[si];
                for (auto [i,j,c]: A)
                    if (c!=0) [[unlikely]]
                        if (i!=j) [[likely]]
                            x[i]-=sol[si][j]*c;
                for (uint32_t i=0;i<A.size().first;i++) 
                    x[i]/=A.ge(i,i);
                r=A*x-b[si];
                mes_r=std::max(mes_r,std::sqrt(r|r));
                sol[si]=x;
            }
        }
        return  sol;
    }

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> GaussZeidelSolver(const Mtype& A, const VMatrix<dtype>& b, uint32_t max_iter=preset_max_iter, dtype max_r=preset_max_r)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();;
        uint32_t iter=0;
        Vector<dtype> x(b.size().first);
        Vector<dtype> r(b.size().first);
        std::vector<std::tuple<uint32_t,uint32_t,dtype>> promise; //promised for whom, required, const
        promise.reserve(b.size().first*4); //4 is a funny constant
        bool first=true;
        while (iter++<max_iter && mes_r>max_r) {
            mes_r=0;
            for (uint32_t si=0;si<b.size().second;si++)
            {
                x=b[si];
                for (auto [i,j,c]: A)
                    if (i<j) x[i]-=sol[si][j]*c;
                    else 
                        if (first) [[unlikely]]
                            if (c!=0) [[unlikely]]
                                if (i!=j) [[likely]]
                                    promise.emplace_back(i,j,c); //will be fast because of smart alloc
                
                if (first) [[unlikely]]
                {
                    std::sort(promise.begin(),promise.end(),std::less<decltype(promise.back())>());
                    first=false;
                }

                for (uint32_t i=0;i<A.size().first;i++)
                    x[i]/=A.ge(i,i);

                for (auto& [i,j,c]: promise)
                {
                    x[i]-=x[j]*(c/A[i,i]);
                }
                r=A*x-b[si];
                mes_r=std::max(mes_r,std::sqrt(r|r));
                sol[si]=x;
            }
        }
        return  sol;
    }


    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> SimpleSolver(const Mtype& A, const VMatrix<dtype>& b, dtype tau, uint32_t max_iter=preset_max_iter, dtype max_r=preset_max_r)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();
        uint32_t iter=0;
        while (iter++<max_iter && mes_r>max_r) {
            mes_r=0;
            for (uint32_t i=0;i<b.size().second;i++)
            {
                Vector<dtype> r=A*sol[i]-b[i];
                sol[i]=sol[i]-r*tau;
                mes_r=std::max(mes_r,std::sqrt(r|r));
            }
        }
        return sol;
    }

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> SteepestGD(const Mtype& A, const VMatrix<dtype>& b, uint32_t max_iter=preset_max_iter, dtype max_r=preset_max_r)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();
        uint32_t iter=0;
        Vector<dtype> r(b.size().first);
        while (iter++<max_iter && mes_r>max_r) {
            mes_r=0;
            for (uint32_t i=0;i<b.size().second;i++)
            {
                r=A*sol[i]-b[i];
                dtype tau=(r|r)/(r|(A*r));
                if (!(tau==tau)) break;
                sol[i]=sol[i]-r*tau;
                mes_r=std::max(mes_r,std::sqrt(r|r));
            }
        }
        return sol;
    }

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> CGD(const Mtype& A, const VMatrix<dtype>& b, uint32_t max_iter=preset_max_iter, dtype max_r=preset_max_r)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();
        uint32_t iter=0;
        VMatrix<dtype> r(b.size()), d(b.size());
        //Vector<dtype> rp(b.size().first);
        dtype alph, beta,rr;
        for (uint32_t i=0;i<b.size().second;i++)
        {
            r[i]=-b[i];
            d[i]=-b[i];
        }
        while (iter++<max_iter && mes_r>max_r) {
            mes_r=0;
            for (uint32_t i=0;i<b.size().second;i++)
            {
                rr=r[i]|r[i];
                //rp=r[i];
                alph=rr/(d[i]|(A*d[i]));
                sol[i]=sol[i]-alph*d[i];
                
                r[i]=A*sol[i]-b[i];
                beta=(r[i]|r[i])/rr;
                
                d[i]=r[i]+beta*d[i];
                //std::cout<<(rp|(r[i]))<<" "<<(d[i]|d[i])<<std::endl;
                mes_r=std::max(mes_r,std::sqrt(r[i]|r[i]));
            }
        }
        return sol;
    }

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    dtype lb_maxIter(const Mtype& A, uint32_t max_iter=preset_max_iter, dtype max_r=preset_max_r)
    {
        Vector<dtype> v(A.size().second);
        for (uint32_t i=0;i<v.size().first;i++)
            v[i]=1;
        Vector<dtype> vold(v.size());
        uint32_t iter=0;
        while ((v-vold).lenght()>max_r && iter++<max_iter) {
            vold=v;
            v=A*v;
            v=v/v.lenght();
        }
        v=A*v;
        return v.lenght();
    }

    template<uint32_t batch,typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> ChebSimpleSolver(const Mtype& A, const VMatrix<dtype>& b, dtype lb_min, dtype lb_max, uint32_t max_iter=preset_max_iter, dtype max_r=preset_max_r)
    {

        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        std::array<dtype, batch> t=LinTrasf<batch>(ChebRootsPermutaton<batch,dtype>(),lb_min,lb_max);
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();
        uint32_t iter=0;
        uint32_t biter=0;
        while (iter++<max_iter && mes_r>max_r) {
            mes_r=0;
            for (uint32_t i=0;i<b.size().second;i++)
            {
                Vector<dtype> r=A*sol[i]-b[i];
                sol[i]=sol[i]-r/t[biter];
                mes_r=std::max(mes_r,std::sqrt(r|r));
            }
            biter=(++biter)&(batch-1);
        }
        return sol;
    }

    template<typename dtype, typename Mtype>
    using StepSig = VMatrix<dtype>(const Mtype&, const VMatrix<dtype>&,const VMatrix<dtype>&);

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> ChebSymAccel(const Mtype& A, const VMatrix<dtype>& b, std::function<StepSig<dtype, Mtype>> step, dtype rho, uint32_t max_iter=preset_max_iter, dtype max_r=preset_max_r)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size()),tmp(b.size());
        sol=step(A,b,tmp);
        dtype mes_r=std::numeric_limits<dtype>::max();
        uint32_t iter=0;
        //std::array<dtype, 3> mu={1,1/rho,2/(rho*rho)-1};
        std::array<dtype, 3> n={1,1/rho,rho-2/rho}; // n[k]=mu[k]/mu[k-1] - want to be near 1
        //n[0] is not used
        while (iter++<max_iter && mes_r>max_r) {
            //                                                     this is bad but faster v
            tmp=std::exchange(sol, (2/(rho*n[2])) * step(A,b,sol) - tmp*(1/(n[2]*n[1])));  
            //mu[0]=mu[1];
            //mu[1]=std::exchange(mu[2], 2*mu[2]/rho-mu[1]);
            n[1]=n[2];
            n[2]=2/rho-1/n[1];
            mes_r=0;
            for (uint32_t i=0;i<b.size().second;i++)
                mes_r=std::max(mes_r,(A*sol[i]-b[i]).lenght());
        }
        //std::cout<<mu[0]<<mu[1]<<mu[2]<<std::endl;
        return sol;
    }

    template<typename dtype, typename Mtype>
    requires IsMatrix<dtype, Mtype>
    VMatrix<dtype> SStepper(const Mtype& A, const VMatrix<dtype>& b, std::function<StepSig<dtype, Mtype>> step, dtype rho, uint32_t max_iter=preset_max_iter, dtype max_r=preset_max_r)
    {
        if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
        VMatrix<dtype> sol(b.size());
        dtype mes_r=std::numeric_limits<dtype>::max();
        uint32_t iter=0;
        while (iter++<max_iter && mes_r>max_r) {
            sol=step(A,b,sol);  
            mes_r=0;
            for (uint32_t i=0;i<b.size().second;i++)
                mes_r=std::max(mes_r,(A*sol[i]-b[i]).lenght());
        }
        //std::cout<<mu[0]<<mu[1]<<mu[2]<<std::endl;
        return sol;
    }
    
}


#endif