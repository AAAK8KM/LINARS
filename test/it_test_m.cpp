#include "imatrix.hpp"
#include "m3diag.hpp"
#include "mdok.hpp"
#include "types.hpp"
#include <iostream>
#include <ostream>
#include <utility>


int main()
{
    Matrix<double> A(std::make_pair(2, 2));
    A.ge(0, 1)=1;
    A.ge(1, 1)=4;
    A.ge(1, 0)=2;
    for (auto [i,j,c] : A)
    {
        std::cout<<i<<" "<<j<<" "<<c<<std::endl;
    }
    std::cout<<std::endl;
    M3diag<double> A2(3);
    A2.ge(0,0)=1;
    A2.ge(0,1)=2;
    A2.ge(1,1)=3;
    A2.ge(1,0)=4;
    A2.ge(2,2)=6;
    A2.ge(2,1)=5;
    for (auto [i,j,c] : A2)
    {
        std::cout<<i<<" "<<j<<" "<<c<<std::endl;
    }
    MDOK<double> A3(A2);
    std::cout<<std::endl;
    for (auto [i,j,c] : A3)
    {
        std::cout<<i<<" "<<j<<" "<<c<<std::endl;
    }

}