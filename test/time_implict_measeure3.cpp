#include "gmres.hpp"
#include "imatrix.hpp"
#include "implictsolver.hpp"
#include "implictstep.hpp"
#include "matrix.hpp"
#include "mcsr.hpp"
#include "mdok.hpp"
#include "mgenerator.hpp"
#include "t2m.hpp"
#include "timer.hpp"
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <random>
#include <utility>
#include <vector>
using namespace LINARS;

std::random_device dev;
std::mt19937 rng(dev());

    constexpr const uint32_t s=1000;
    constexpr const uint32_t size=s*s, max_test_it=50, num_t=100;

std::fstream file;

template <uint32_t N>
void some(const MCSR<long double>& A, const VMatrix<long double>& b, const VMatrix<long double>& exp)
{
      
   
    if (!file.is_open()) {
        std::cerr << "Error creating file"  << std::endl;
        return;
    }
    //exit(0);
    std::cout<<"start"<<N<<std::endl;
    auto gs = [](const MCSR<long double>& A, const VMatrix<long double>& b,const VMatrix<long double>& prev)->VMatrix<long double>{
        return GMRES<N,long double,MCSR<long double>>(A,b,prev);
    };
    auto [res,duration] = ms_timer(SStepper<long double,MCSR<long double>>, A, b, std::function<StepSig<long double, MCSR<long double>>>(gs), num_t, 0.);
    auto r1=sqrt((res[0]-exp[0])|(res[0]-exp[0]));
    auto r2=sqrt((A*res[0]-b[0])|(A*res[0]-b[0]));
    file<<N<<","<<r1<<","<<norm2(res[0]-exp[0])<<","<<r2<<","<<duration<<std::endl;
        //std::cout<<res[0]<<std::endl;
    std::cout<<"end";
    //file.close();
}


int main()
{
     file.open("gmresm.csv",  std::ios_base::out);
        //constexpr const uint32_t test_it_step=max_test_it/max_test_num;
    VMatrix<long double> b(size,1), exp(size,1);
    std::uniform_real_distribution<long double> dist(0.1,0.9);
    std::uniform_real_distribution<long double> dist2(-1,1);
    /*for (std::size_t i=0;i<size;i++)
        for (std::size_t j=i;j<size;j++)
            if (dist2(rng)>0.5)
            {
                if (i%2)
                    A.set(i, j,dist(rng)*10);
                else
                    A.set(j,i,dist(rng)*10);
            }
    for (std::size_t i=0;i<size;i++)
    {
        A.set(i, i, dist(rng));
        exp[0][i]=dist(rng)*10;
    }*/
    for (std::size_t i=0;i<size;i++)
    {
        exp[0][i]=sin((i+0.)/s+i%s);
    }
    MCSR<long double> A=PuassonTask0<long double,MCSR<long double>>(s,s);
    b[0]=A*exp[0];
    //std::cout<<A<<b<<std::endl;//<<exp<<std::endl;
    some<1>(A, b, exp);
    some<2>(A, b, exp);
    some<3>(A, b, exp);
    some<4>(A, b, exp);
    some<5>(A, b, exp);
    some<6>(A, b, exp);
    some<7>(A, b, exp);
    some<8>(A, b, exp);
    some<9>(A, b, exp);
    some<10>(A, b, exp);
    
    return 0;
    
}