#include "imatrix.hpp"
#include "matrix.hpp"
#include "mcsr.hpp"
#include "mdok.hpp"
#include "mgenerator.hpp"
#include "implictstep.hpp"
#include "implictsolver.hpp"
#include "qrdec.hpp"
#include "qrsolver.hpp"
#include "trmatrix.hpp"
#include "lvlmatrix.hpp"
#include "holetski.hpp"
#include "gmres.hpp"
#include <fstream>
#include "t2m.hpp"
#include "mvector.hpp"
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <ostream>
#include <utility>

using namespace LINARS;


/*int main()
{
    MDOK<double> A(3,3);
    A[0,1]=1;
    A[0,2]=2;
    A[1,0]=1;
    A[1,2]=2;
    A[2,0]=2;
    A[2,2]=1;
    //std::cout<<A;
    Vector<double> b(3),prev(3);
    b[2]=1;
    //std::cout<<b; 
    auto D=GMRESdata<2,double>(b.size().first);
    GMRES_R<2,double,MDOK<double>>(A,b,prev,D);
    std::cout<<D.basis<<std::endl;
    std::cout<<D.R<<D.r;
}*/

/*int main()
{
    MDOK<double> A(3,3);
    A[0,0]=10;
    A[0,1]=1;
    A[1,1]=5;
    A[2,0]=1;
    A[2,2]=2;
    std::cout<<lvlmatrix<double>(A);
}*/

int main()
{
    MDOK<long double> A(3,3);
    A[0,0]=10;
    A[0,1]=3;
    A[0,2]=6;
    A[1,0]=3;
    A[1,1]=5;
    A[1,2]=1;
    A[2,0]=6;
    A[2,1]=1;
    A[2,2]=8;
    Vector<long double> b(3);
    b[0]=1;
    b[1]=2;
    b[2]=1;
    
    SteepestGD<long double,MDOK<long double>>(A, b);
    CGD<long double,MDOK<long double>>(A, b);
    SimpleSolver<long double,MDOK<long double>>(A, b, 0.11029174063629871);
}