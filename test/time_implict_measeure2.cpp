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
    constexpr const uint32_t s=500;
    constexpr const uint32_t size=s*s, max_test_it=50, max_test_num=50;
    constexpr const uint32_t test_it_step=max_test_it/max_test_num;
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
    std::fstream file;
    /*file.open("jacobi.csv",  std::ios_base::out);
    if (!file.is_open()) {
        std::cerr << "Error creating file"  << std::endl;
        return 1;
    }
    std::cout<<"start"<<std::endl;
    for (uint64_t num_t=test_it_step;num_t<=max_test_it;num_t+=test_it_step)
    {
        //JakobiSolver - 79  GaussZeidelSolver - 40  SimpleSolver
        auto [res,duration] = ms_timer(JakobiSolver<long double,MCSR<long double>>,A, b,num_t,0.);
        auto r1=sqrt((res[0]-exp[0])|(res[0]-exp[0]));
        auto r2=sqrt((A*res[0]-b[0])|(A*res[0]-b[0]));
        file<<num_t<<","<<r1<<","<<norm2(res[0]-exp[0])<<","<<r2<<","<<duration<<std::endl;
        //std::cout<<res[0]<<std::endl;
        std::cout<<"work "<<(num_t*100.)/max_test_it<<"%"<<std::endl;
    }
    file.close();*/
    file.open("gauszeidel.csv",  std::ios_base::out);
    if (!file.is_open()) {
        std::cerr << "Error creating file"  << std::endl;
        return 1;
    }
    std::cout<<"start"<<std::endl;
    for (uint64_t num_t=test_it_step;num_t<=max_test_it;num_t+=test_it_step)
    {
        //JakobiSolver - 79  GaussZeidelSolver - 40  SimpleSolver
        auto [res,duration] = ms_timer(GaussZeidelSolver<long double,MCSR<long double>>,A, b,num_t,0.);
        auto r1=sqrt((res[0]-exp[0])|(res[0]-exp[0]));
        auto r2=sqrt((A*res[0]-b[0])|(A*res[0]-b[0]));
        file<<num_t<<","<<r1<<","<<norm2(res[0]-exp[0])<<","<<r2<<","<<duration<<std::endl;
        //std::cout<<res[0]<<std::endl;
        std::cout<<"work "<<(num_t*100.)/max_test_it<<"%"<<std::endl;
    }
    file.close();
    long double lb_max=lb_maxIter<long double>(A,100), lb_min=18./((s+1)*(s+1));
    std::cout<<lb_max<<std::endl;

    /*file.open("simple.csv",  std::ios_base::out);
    if (!file.is_open()) {
        std::cerr << "Error creating file"  << std::endl;
        return 1;
    }
    std::cout<<"start"<<std::endl;
    for (uint64_t num_t=test_it_step;num_t<=max_test_it;num_t+=test_it_step)
    {
        //JakobiSolver - 79  GaussZeidelSolver - 40  SimpleSolver
        auto [res,duration] = ms_timer(SimpleSolver<long double,MCSR<long double>>,A, b,2/(lb_max+lb_min),num_t,0.);
        auto r1=sqrt((res[0]-exp[0])|(res[0]-exp[0]));
        auto r2=sqrt((A*res[0]-b[0])|(A*res[0]-b[0]));
        file<<num_t<<","<<r1<<","<<norm2(res[0]-exp[0])<<","<<r2<<","<<duration<<std::endl;
        //std::cout<<res[0]<<std::endl;
        std::cout<<"work "<<(num_t*100.)/max_test_it<<"%"<<std::endl;
    }
    file.close();*/
    file.open("SteeperGD.csv",  std::ios_base::out);
    if (!file.is_open()) {
        std::cerr << "Error creating file"  << std::endl;
        return 1;
    }
    std::cout<<"start"<<std::endl;
    for (uint64_t num_t=test_it_step;num_t<=max_test_it;num_t+=test_it_step)
    {
        //JakobiSolver - 79  GaussZeidelSolver - 40  SimpleSolver
        auto [res,duration] = ms_timer(SteepestGD<long double,MCSR<long double>>,A, b,num_t,0.);
        auto r1=sqrt((res[0]-exp[0])|(res[0]-exp[0]));
        auto r2=sqrt((A*res[0]-b[0])|(A*res[0]-b[0]));
        file<<num_t<<","<<r1<<","<<norm2(res[0]-exp[0])<<","<<r2<<","<<duration<<std::endl;
        //std::cout<<res[0]<<std::endl;
        std::cout<<"work "<<(num_t*100.)/max_test_it<<"%"<<std::endl;
    }
    file.close();
    file.open("CGD.csv",  std::ios_base::out);
    if (!file.is_open()) {
        std::cerr << "Error creating file"  << std::endl;
        return 1;
    }
    std::cout<<"start"<<std::endl;
    for (uint64_t num_t=test_it_step;num_t<=max_test_it;num_t+=test_it_step)
    {
        //JakobiSolver - 79  GaussZeidelSolver - 40  SimpleSolver
        auto [res,duration] = ms_timer(CGD<long double,MCSR<long double>>,A, b,num_t,0.);
        auto r1=sqrt((res[0]-exp[0])|(res[0]-exp[0]));
        auto r2=sqrt((A*res[0]-b[0])|(A*res[0]-b[0]));
        file<<num_t<<","<<r1<<","<<norm2(res[0]-exp[0])<<","<<r2<<","<<duration<<std::endl;
        //std::cout<<res[0]<<std::endl;
        std::cout<<"work "<<(num_t*100.)/max_test_it<<"%"<<std::endl;
    }
    file.close();
    /*file.open("simplecheb.csv",  std::ios_base::out);
    if (!file.is_open()) {
        std::cerr << "Error creating file"  << std::endl;
        return 1;
    }
    //exit(0);
    std::cout<<"start"<<std::endl;
    for (uint64_t num_t=test_it_step;num_t<=max_test_it;num_t+=test_it_step)
    {
        //JakobiSolver - 79  GaussZeidelSolver - 40  SimpleSolver
        auto [res,duration] = ms_timer(ChebSimpleSolver<8,long double,MCSR<long double>>,A, b,lb_min,lb_max,num_t,0.);
        auto r1=sqrt((res[0]-exp[0])|(res[0]-exp[0]));
        auto r2=sqrt((A*res[0]-b[0])|(A*res[0]-b[0]));
        file<<num_t<<","<<r1<<","<<norm2(res[0]-exp[0])<<","<<r2<<","<<duration<<std::endl;
        //std::cout<<res[0]<<std::endl;
        std::cout<<"work "<<(num_t*100.)/max_test_it<<"%"<<std::endl;
    }
    file.close();*/
    auto ls = [lb_max,p=SPrep<long double,MCSR<long double>>(A)](const MCSR<long double>& A, const VMatrix<long double>& b,const VMatrix<long double>& prev)->VMatrix<long double>{
        return SSORStep<long double,MCSR<long double>>(A,b,prev,p,1+1/(4*lb_max*lb_max));
    };
    auto l = [lb_max,p=SPrep<long double,MCSR<long double>>(A)](const MCSR<long double>& A, const VMatrix<long double>& b,const VMatrix<long double>& prev)->VMatrix<long double>{
        return SORStep<long double,MCSR<long double>>(A,b,prev,p,1+1/(4*lb_max*lb_max));
    };

    file.open("ssorchb.csv",  std::ios_base::out);
    if (!file.is_open()) {
        std::cerr << "Error creating file"  << std::endl;
        return 1;
    }
    //exit(0);
    std::cout<<"start"<<std::endl;
    for (uint64_t num_t=test_it_step;num_t<=max_test_it;num_t+=test_it_step)
    {
        auto [res,duration] = ms_timer(ChebSymAccel<long double,MCSR<long double>>, A, b, std::function<StepSig<long double, MCSR<long double>>>(ls), 1/lb_max, num_t, 0.);
        auto r1=sqrt((res[0]-exp[0])|(res[0]-exp[0]));
        auto r2=sqrt((A*res[0]-b[0])|(A*res[0]-b[0]));
        file<<num_t<<","<<r1<<","<<norm2(res[0]-exp[0])<<","<<r2<<","<<duration<<std::endl;
        //std::cout<<res[0]<<std::endl;
        std::cout<<"work "<<(num_t*100.)/max_test_it<<"%"<<std::endl;
    }
    file.close();
    /*file.open("ssor.csv",  std::ios_base::out);
    if (!file.is_open()) {
        std::cerr << "Error creating file"  << std::endl;
        return 1;
    }
    //exit(0);
    std::cout<<"start"<<std::endl;
    for (uint64_t num_t=test_it_step;num_t<=max_test_it;num_t+=test_it_step)
    {
        auto [res,duration] = ms_timer(SStepper<long double,MCSR<long double>>, A, b, std::function<StepSig<long double, MCSR<long double>>>(ls), 1/lb_max, num_t, 0.);
        auto r1=sqrt((res[0]-exp[0])|(res[0]-exp[0]));
        auto r2=sqrt((A*res[0]-b[0])|(A*res[0]-b[0]));
        file<<num_t<<","<<r1<<","<<norm2(res[0]-exp[0])<<","<<r2<<","<<duration<<std::endl;
        //std::cout<<res[0]<<std::endl;
        std::cout<<"work "<<(num_t*100.)/max_test_it<<"%"<<std::endl;
    }
    file.close();*/
    file.open("sor.csv",  std::ios_base::out);
    if (!file.is_open()) {
        std::cerr << "Error creating file"  << std::endl;
        return 1;
    }
    //exit(0);
    std::cout<<"start"<<std::endl;
    for (uint64_t num_t=test_it_step;num_t<=max_test_it;num_t+=test_it_step)
    {
        auto [res,duration] = ms_timer(SStepper<long double,MCSR<long double>>, A, b, std::function<StepSig<long double, MCSR<long double>>>(l), 1/lb_max, num_t, 0.);
        auto r1=sqrt((res[0]-exp[0])|(res[0]-exp[0]));
        auto r2=sqrt((A*res[0]-b[0])|(A*res[0]-b[0]));
        file<<num_t<<","<<r1<<","<<norm2(res[0]-exp[0])<<","<<r2<<","<<duration<<std::endl;
        //std::cout<<res[0]<<std::endl;
        std::cout<<"work "<<(num_t*100.)/max_test_it<<"%"<<std::endl;
    }
    file.close();
}