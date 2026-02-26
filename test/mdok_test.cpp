#include "mdok.hpp"
#include <cstdint>
#include <random>
#include <utility>
#include <gtest/gtest.h>

TEST(DOK_Matrix_Test, manual)
{
    Matrix<double> A2(std::make_pair(5, 5));
    A2.ge(0,0)=1;
    A2.ge(0,2)=2;
    A2.ge(1,3)=3;
    A2.ge(1,1)=4;
    A2.ge(2,1)=5;
    A2.ge(4,4)=6;
    EXPECT_EQ(A2,  MDOK<double> (A2));
    EXPECT_EQ(A2,  Matrix<double>(MDOK<double> (A2)));
}

class Mtest : public ::testing::TestWithParam<int> {};

std::random_device dev;
std::mt19937 rng(dev());

TEST_P(Mtest, random)
{
    Matrix<double> A2(std::make_pair(100, 100));
    std::uniform_real_distribution<float> dist1(0,1);
    std::uniform_real_distribution<float> dist2(-1000,1000);
    float c=dist1(rng);
    
    for (uint32_t i=0;i<10000;i++)
        if (c<dist1(rng))
            A2.ge(i/100, i%100)=dist2(rng);

    EXPECT_EQ(A2,   MDOK<double> (A2));
}

INSTANTIATE_TEST_SUITE_P(
    DOK_Matrix_Test, 
    Mtest,   
    ::testing::Range(0, 10) 
);