#include "gmres.hpp"
#include "holetski.hpp"
#include "imatrix.hpp"
#include "implictsolver.hpp"
#include "implictstep.hpp"
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

int main()
{
    constexpr const uint32_t s=50;
    constexpr const uint32_t size=s*s, max_test_it=50, max_test_num=50;
    constexpr const uint32_t test_it_step=max_test_it/max_test_num;
    VMatrix<long double> b(size,1), exp(size,1);
    std::uniform_real_distribution<long double> dist(0.1,0.9);
    std::uniform_real_distribution<long double> dist2(-1,1);
    for (std::size_t i=0;i<size;i++)
    {
        exp[0][i]=sin((i+0.)/s+i%s);
    }
    MCSR<long double> A=PuassonTask0<long double,MCSR<long double>>(s,s);
    b[0]=A*exp[0];
    //std::cout<<A<<b<<std::endl;//<<exp<<std::endl;
    std::fstream file;

    file.open("CGD.csv",  std::ios_base::out);
    if (!file.is_open()) {
        std::cerr << "Error creating file"  << std::endl;
        return 1;
    }
    std::cout<<"start"<<std::endl;
    for (uint64_t num_t=test_it_step;num_t<=max_test_it;num_t+=test_it_step)
    {
        //JakobiSolver - 79  GaussZeidelSolver - 40  SimpleSolver
        auto [res,duration] = ms_timer<20>(CGD<long double,MCSR<long double>>,A, b,num_t,0.);
        auto r1=sqrt((res[0]-exp[0])|(res[0]-exp[0]));
        auto r2=sqrt((A*res[0]-b[0])|(A*res[0]-b[0]));
        file<<num_t<<","<<r1<<","<<norm2(res[0]-exp[0])<<","<<r2<<","<<duration<<std::endl;
        //std::cout<<res[0]<<std::endl;
        std::cout<<"work "<<(num_t*100.)/max_test_it<<"%"<<std::endl;
    }
    file.close();

    file.open("PCCGD.csv",  std::ios_base::out);
    if (!file.is_open()) {
        std::cerr << "Error creating file"  << std::endl;
        return 1;
    }
    std::cout<<"start"<<std::endl;
    auto L=holetski<long double>(A);
    for (uint64_t num_t=test_it_step;num_t<=max_test_it;num_t+=test_it_step)
    {
        //JakobiSolver - 79  GaussZeidelSolver - 40  SimpleSolver
        auto [res,duration] = ms_timer<20>(PC_CGD<long double,MCSR<long double>>,A, b, L,num_t,0.);
        auto r1=sqrt((res[0]-exp[0])|(res[0]-exp[0]));
        auto r2=sqrt((A*res[0]-b[0])|(A*res[0]-b[0]));
        file<<num_t<<","<<r1<<","<<norm2(res[0]-exp[0])<<","<<r2<<","<<duration<<std::endl;
        //std::cout<<res[0]<<std::endl;
        std::cout<<"work "<<(num_t*100.)/max_test_it<<"%"<<std::endl;
    }
    file.close();

}