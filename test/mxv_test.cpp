#include "matrix.hpp"
#include "types.hpp"
#include <cstdint>
#include <random>
#include <utility>
#include <gtest/gtest.h>



std::random_device dev;
std::mt19937 rng(dev());

class MtestDOK : public ::testing::TestWithParam<int> {};

TEST_P(MtestDOK, random)
{
    Matrix<double> A2(std::make_pair(1000, 1000));
    std::uniform_real_distribution<float> dist1(0,1);
    std::uniform_real_distribution<float> dist2(-1000,1000);
    float c=dist1(rng);
    
    Vector<double> v(1000);

    for (uint32_t i=0;i<1000000;i++)
        if (c<dist1(rng))
            A2.ge(i/1000, i%1000)=dist2(rng);
    for (uint32_t i=0;i<1000;i++)
        v[i]=dist2(rng);

    MDOK<double> A3(A2);

    EXPECT_EQ(A3*v,   A2*v);
}

INSTANTIATE_TEST_SUITE_P(
    DOK_Vector_Test, 
    MtestDOK,   
    ::testing::Range(0, 10) 
);


class MtestCSR : public ::testing::TestWithParam<int> {};

TEST_P(MtestCSR, random)
{
    Matrix<double> A2(std::make_pair(1000, 1000));
    std::uniform_real_distribution<float> dist1(0,1);
    std::uniform_real_distribution<float> dist2(-1000,1000);
    float c=dist1(rng);
    
    Vector<double> v(1000);

    for (uint32_t i=0;i<1000000;i++)
        if (c<dist1(rng))
            A2.ge(i/1000, i%1000)=dist2(rng);
    for (uint32_t i=0;i<1000;i++)
        v[i]=dist2(rng);

    MCSR<double> A3(A2);

    EXPECT_EQ(A3*v,   A2*v);
}

INSTANTIATE_TEST_SUITE_P(
    MCSR_Vector_Test, 
    MtestCSR,   
    ::testing::Range(0, 10) 
);