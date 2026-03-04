#include "imatrix.hpp"
#include <cstdint>
#include <random>
#include <utility>
#include <gtest/gtest.h>

using namespace LINARS;



std::random_device dev;
std::mt19937 rng(dev());


class VtestS : public ::testing::TestWithParam<int> {};
TEST_P(VtestS, random)
{
    Vector<double> V1(std::make_pair(100, 1)),V2(V1),VS(V1);
    std::uniform_real_distribution<double> dist2(-1000,1000);
    
    for (uint32_t i=0;i<V1.size().first;i++)
    {
        V1[i]=dist2(rng);
        V2[i]=dist2(rng);
        VS[i]=V1[i]+V2[i];
    }
       

    EXPECT_EQ(V1+V2,  VS);
}
INSTANTIATE_TEST_SUITE_P(
    Vector_Sum_Test, 
    VtestS,   
    ::testing::Range(0, 10) 
);

class VtestMVxD : public ::testing::TestWithParam<int> {};
TEST_P(VtestMVxD, random)
{
    Vector<double> V1(std::make_pair(100, 1)),V2(V1);//,VS(V1);
    std::uniform_real_distribution<double> dist2(-10,10);
    double c=dist2(rng);
    for (uint32_t i=0;i<V1.size().first;i++)
    {
        V1[i]=dist2(rng);
        V2[i]=V1[i]*c;
        //VS[i]=V1[i]+V2[i];
    }
       

    EXPECT_EQ(V1*c,  V2);
    EXPECT_EQ(c*V1,  V2);
}
INSTANTIATE_TEST_SUITE_P(
    Vector_VectorxValue_Test, 
    VtestMVxD,   
    ::testing::Range(0, 10) 
);

class VtestMVdD : public ::testing::TestWithParam<int> {};
TEST_P(VtestMVdD, random)
{
    Vector<double> V1(std::make_pair(100, 1)),V2(V1);//,VS(V1);
    std::uniform_real_distribution<double> dist2(-10,10);
    double c=dist2(rng);
    for (uint32_t i=0;i<V1.size().first;i++)
    {
        V1[i]=dist2(rng);
        V2[i]=V1[i]/c;
        //VS[i]=V1[i]+V2[i];
    }
       

    EXPECT_EQ(V1/c,  V2);
}
INSTANTIATE_TEST_SUITE_P(
    Vector_VectordValue_Test, 
    VtestMVdD,   
    ::testing::Range(0, 10) 
);

class VtestMVxV1 : public ::testing::TestWithParam<int> {};
TEST_P(VtestMVxV1, random)
{
    Vector<double> V1(std::make_pair(100, 1)),V2(std::make_pair(1, 100));//,VS(V1);
    std::uniform_real_distribution<double> dist2(-10,10);
    double c=0;
    for (uint32_t i=0;i<V1.size().first;i++)
    {
        V1[i]=dist2(rng);
        V2[i]=dist2(rng);
        c+=V1[i]*V2[i];
        //VS[i]=V1[i]+V2[i];
    }
    EXPECT_EQ(V2*V1,  c);
}
INSTANTIATE_TEST_SUITE_P(
    Vector_VectorxVector_Test, 
    VtestMVxV1,   
    ::testing::Range(0, 10) 
);

class VtestMVxV2 : public ::testing::TestWithParam<int> {};
TEST_P(VtestMVxV2, random)
{
    Vector<double> V1(std::make_pair(100, 1)),V2(V1);//,VS(V1);
    std::uniform_real_distribution<double> dist2(-10,10);
    double c=0;
    for (uint32_t i=0;i<V1.size().first;i++)
    {
        V1[i]=dist2(rng);
        V2[i]=dist2(rng);
        c+=V1[i]*V2[i];
        //VS[i]=V1[i]+V2[i];
    }
    EXPECT_EQ(V1.transposed_shared()*V2,  c);
    EXPECT_EQ(V2.transposed()*V1,  c);
}
INSTANTIATE_TEST_SUITE_P(
    Vector_VectorxVector2_Test, 
    VtestMVxV2,   
    ::testing::Range(0, 10) 
);

