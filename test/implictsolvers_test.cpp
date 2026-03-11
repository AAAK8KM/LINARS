#include "imatrix.hpp"
#include "implictsolver.hpp"
#include "t2m.hpp"
#include <iostream>
#include <gtest/gtest.h>
#include <ostream>
#include <random>
#include <utility>
using namespace LINARS;

TEST(ImplictSolversTest,Manual)
{
    Matrix<double> A(3,3);
    A.ge(0, 0)=10;
    A.ge(0, 1)=2;
    A.ge(0, 2)=3;
    A.ge(1, 0)=4;
    A.ge(1, 1)=20;
    A.ge(1, 2)=6;
    A.ge(2, 0)=7;
    A.ge(2, 1)=8;
    A.ge(2, 2)=30;
    Vector<double> exp(3);
    exp[0]=1;
    exp[1]=2;
    exp[2]=3;
    Vector<double> b=A*exp;

    auto res1=GaussZeidelSolver<double>(A,b);
    auto res2=JakobiSolver<double>(A,b, 100000);
    auto res3=SimpleSolver<double>(A,b,0.01);

    //std::cout<<res1<<std::endl;

    std::cout<<"GaussZeidelSolver"<<std::endl;
    std::cout<<res1[0]<<std::endl;
    EXPECT_TRUE(((res1[0]-exp)|(res1[0]-exp))<1e-6);
    std::cout<<"JakobiSolver"<<std::endl;
    std::cout<<res2[0]<<std::endl;
    EXPECT_TRUE(((res2[0]-exp)|(res2[0]-exp))<1e-6);
    std::cout<<"SimpleSolver"<<std::endl;
    std::cout<<res3[0]<<std::endl;
    EXPECT_TRUE(((res3[0]-exp)|(res3[0]-exp))<1e-6);
}

/*
TEST(QRSolverTest,Manual2)
{
    Matrix<double> A(3,3);
    A.ge(0, 0)=1;
    A.ge(0, 1)=2;
    A.ge(0, 2)=3;
    A.ge(1, 0)=4;
    A.ge(1, 1)=5;
    A.ge(1, 2)=6;
    A.ge(2, 0)=7;
    A.ge(2, 1)=8;
    A.ge(2, 2)=1;
    VMatrix<double> exp(std::make_pair(3,2));
    exp.ge(0, 0)=1;
    exp.ge(1, 0)=2;
    exp.ge(2, 0)=3;
    exp.ge(0, 1)=4;
    exp.ge(1, 1)=5;
    exp.ge(2, 1)=6;
    auto b=A*exp;

    auto res=QRSolver<double>(A,b,QRdecompositionH);

    std::cout<<typeid(res).name()<<" "<<res.size().first<<" " <<res.size().second<<std::endl;
    std::cout<<res<<std::endl;

    EXPECT_EQ(exp, res);
}

std::random_device dev;
std::mt19937 rng(dev());

class QRtestHSolve : public ::testing::TestWithParam<int> {};
TEST_P(QRtestHSolve, random)
{
    Matrix<double> A(100,100);
    Matrix<double> B(A);
    std::uniform_real_distribution<double> dist2(-10,10);
    for (auto [i,j,c]: A)
    {
        A.ge(i, j)=dist2(rng);
        B.ge(i, j)=dist2(rng);
    }

    EXPECT_EQ(B,  QRSolver<double>(A,A*B,QRdecompositionH));
}
INSTANTIATE_TEST_SUITE_P(
    QRHaushold_solve_Test, 
    QRtestHSolve,   
    ::testing::Range(0, 10) 
);*/