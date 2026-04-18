#include "imatrix.hpp"
#include "mcsr.hpp"
#include "mgenerator.hpp"
#include "implictstep.hpp"
#include "implictsolver.hpp"
#include "qrdec.hpp"
#include "qrsolver.hpp"
#include <fstream>
#include "t2m.hpp"
#include "mvector.hpp"
#include <cstdint>
#include <functional>
#include <iostream>
#include <ostream>
#include <utility>

using namespace LINARS;

/*int main()
{
    Matrix<double>A(4,4);
    A.ge(0, 1)=8;
    A.ge(1, 0)=8;
    A.ge(1, 3)=2;
    A.ge(2, 2)=4;
    A.ge(3, 1)=1;
    A.ge(3, 2)=2;
    std::cout<<A<<std::endl;
    MCSR<double> M(A);
    M.flush();
}*/

/*int main()
{
    Matrix<double>A(4,4);
    A.ge(0, 0)=90;
    A.ge(1, 1)=65;
    A.ge(2, 2)=77;
    A.ge(3, 3)=50;
    A.ge(0, 1)=4;
    A.ge(1, 2)=8;
    A.ge(1, 3)=4;
    A.ge(2, 3)=5;
    A.ge(3, 0)=5;
    A.ge(3, 1)=4;
    Matrix<double>B(4,4);
    B.ge(0, 0)=0.91;
    B.ge(1, 1)=0.935;
    B.ge(2, 2)=0.923;
    B.ge(3, 3)=0.95;
    B.ge(0, 1)=-4./8000;
    B.ge(1, 2)=-8./8000;
    B.ge(1, 3)=-4./8000;
    B.ge(2, 3)=-5./8000;
    B.ge(3, 0)=-5./8000;
    B.ge(3, 1)=-4./8000;
    Vector<double>b(4);
    b[0]=6;
    b[1]=6;
    b[2]=7;
    b[3]=4;
    std::cout<<A<<std::endl<<b<<std::endl;
    Vector<double> c=QRSolver<double>(A, b, QRdecompositionH<double>);
    std::cout<<c<<std::endl; 
    std::cout<<c.lenght()<<std::endl;  
    std::cout<<B<<std::endl; 
    std::cout<<lb_maxIter<double>(B)<<std::endl;   
    
}*/
/*
template<typename dtype, typename Mtype>
requires IsMatrix<dtype, Mtype>
std::pair<uint32_t, bool> SimpleSolverT(const Mtype& A, const Vector<dtype>& b, dtype tau, dtype max_r=1e-12)
{
    if (A.size().first!=b.size().first) throw std::runtime_error("Invalid linar system");
    VMatrix<dtype> sol(b.size());
    dtype mes_r=std::numeric_limits<dtype>::max();
    uint32_t iter=0;
    while (mes_r>max_r && (mes_r<1e2 || iter==0)) {
        mes_r=0;
        for (uint32_t i=0;i<b.size().second;i++)
        {
            Vector<dtype> r=A*sol[i]-b;
            sol[i]=sol[i]-r*tau;
            mes_r=std::max(mes_r,std::sqrt(r|r));
        }
        iter++;
    }
    return std::make_pair(iter,mes_r<1e16);
}


int main(){
    Matrix<long double>A(4,4);
    A.ge(0, 0)=79;//80;
    A.ge(1, 1)=80;//80.;
    A.ge(2, 2)=56;//80.;
    A.ge(3, 3)=57;//80.;
    A.ge(0, 1)=2;//-80.;
    A.ge(0, 2)=6;//-80.;
    A.ge(1, 0)=2;//-80.;
    A.ge(2, 1)=3;//-80.;
    A.ge(3, 0)=8;//-80.;
    Vector<long double>b(4);
    b[0]=1;
    b[1]=1;
    b[2]=1;
    b[3]=1;
    std::cout<<A<<b<<std::endl;
    std::fstream file("tau4.csv", std::ios_base::out);
    bool flag=1;
    long double tau=1;
    while (flag && tau<=1000) {
        auto [N, res]=SimpleSolverT(A,b,tau);
        file<<tau<<" "<<N<<std::endl;
        flag=res;
        tau+=(1000.-1)/1000;
    }
    
}*/


/*int main()
{
    Matrix<double> A(3,3);
    A[0, 0]=10;
    A[1, 1]=10;
    A[2, 2]=10;
    A[2, 0]=2;
    A[0, 2]=1;
    VMatrix<double> b(3,1);
    b[0, 0]=1;
    b[1, 0]=2;
    b[2, 0]=3;
    VMatrix<double> x(3,1);

    for (auto& pr: SSORPrep<double,Matrix<double>>(A))
    {
        std::cout<<"a\n";
        for (auto [i,j,c]: pr)
            std::cout<<i<<" "<<j<<" "<<c<<std::endl;
        std::cout<<std::endl;
    }
    auto l = [p=SSORPrep<double,Matrix<double>>(A)](const Matrix<double>& A, const VMatrix<double>& b,const VMatrix<double>& prev)->VMatrix<double>{
        return SSORStep<double,Matrix<double>>(A,b,prev,p);
    };
    x=ChebSymAccel(A,b,std::function<StepSig<double, Matrix<double>>>(l),0.9);
    std::cout<<x<<std::endl;

}*/

int main()
{
    //MCSR<double> M = PuassonTask0<double,MCSR<double>>(5,5);
    std::cout<<"strat\n";
    Vector<double> M(5);
    M=M+M;
    std::cout<<M;
}