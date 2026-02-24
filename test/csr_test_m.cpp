#include "imatrix.hpp"
#include "matrix.hpp"
#include <cstdint>
#include <iostream>
#include <ostream>
#include <utility>


int main()
{
    Matrix<double> A2(std::make_pair(5, 5));
    A2.ge(0,0)=1;
    A2.ge(0,2)=2;
    A2.ge(1,3)=3;
    A2.ge(1,1)=4;
    A2.ge(2,1)=5;
    A2.ge(4,4)=6;
    std::cout<<"ok0"<<std::endl;
    MCSR<double> A3(A2);
    for (uint32_t i=0;i<5;i++)
    {
        for (uint32_t j=0;j<5;j++)
            std::cout<<A3.gev(i, j)<<" ";
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    Matrix<double> A4(A3);
    for (uint32_t i=0;i<5;i++)
    {
        for (uint32_t j=0;j<5;j++)
            std::cout<<A4.gev(i, j)<<" ";
        std::cout<<std::endl;
    }
}