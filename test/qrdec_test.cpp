#include "imatrix.hpp"
#include "matrix.hpp"
#include "t2m.hpp"
#include "qrdec.hpp"
#include <iostream>
#include <ostream>
#include <gtest/gtest.h>

using namespace LINARS;

TEST(QRDecomposition,hand_test)
{
    Matrix<double> E(3,3);
    E.ge(0, 0)=1;
    E.ge(1, 1)=1;
    E.ge(2, 2)=1;
    Matrix<double> A(3,3);
    A.ge(0, 0)=3;
    A.ge(0, 1)=10;
    A.ge(0, 2)=17;
    A.ge(1, 0)=4;
    A.ge(1, 1)=5;
    A.ge(1, 2)=6;
    A.ge(2, 2)=5;
    std::cout<<A<<std::endl;
    auto [Q,R] = QRdecompositionH<double>(A);
    Q.reverse();
    std::cout<<R<<std::endl<<Q(E)<<std::endl<<Q.opt(A)<<std::endl;
    EXPECT_EQ(Q(A), R);
    EXPECT_EQ(Q.opt(A), R);
    Q.reverse();
    std::cout<<Q(R);
    EXPECT_EQ(Q(R), A);
    EXPECT_EQ(Q.opt(R), A);

}

TEST(QRDecomposition,hand_test2)
{
    Matrix<double> E(3,3);
    E.ge(0, 0)=1;
    E.ge(1, 1)=1;
    E.ge(2, 2)=1;
    Matrix<double> A(3,3);
    A.ge(0, 0)=1;
    A.ge(0, 1)=2;
    A.ge(0, 2)=3;
    A.ge(1, 0)=4;
    A.ge(1, 1)=5;
    A.ge(1, 2)=6;
    A.ge(2, 0)=7;
    A.ge(2, 1)=8;
    A.ge(2, 2)=9;
    std::cout<<A<<std::endl;
    auto [Q,R] = QRdecompositionH<double>(A);
    std::cout<<"check res\n";
    std::cout<<R<<std::endl<<Q(E)<<std::endl<<Q.opt(A)<<std::endl;
    Q.reverse();
    EXPECT_EQ(Q(A), R);
    EXPECT_EQ(Q.opt(A), R);
    Q.reverse();
    std::cout<<Q(R);
    EXPECT_EQ(Q(R), A);
    EXPECT_EQ(Q.opt(R), A);

}