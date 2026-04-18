#include "matrixes.hpp"
#include "t2m.hpp"
#include "tomas.hpp"
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <ostream>
#include <random>
#include <stdexcept>

using namespace LINARS;

std::random_device dev;
std::mt19937 rng(dev());


template<typename dtype>
void check(std::size_t N, std::size_t k)
{
    std::uniform_real_distribution<dtype> dist(1,10);
    M3diag<dtype> A(N);
    Matrix<dtype> b(N,1), res(N,1);
    for (std::size_t i=0;i<k;i++)
    {
        for (std::size_t j=0;j<N-1;j++)
        {
            A.ge(j+1, j)=dist(rng);
            A.ge(j, j+1)=dist(rng);
        }
        dtype d;
        for (std::size_t j=0;j<N;j++)
        {
            d=0;
            if (j<N-1)
                d+=std::max(std::abs(A.ge(j+1, j)),std::abs(A.ge(j, j+1)));
            if (j>0)
                d+=std::max(std::abs(A.ge(j-1, j)),std::abs(A.ge(j, j-1)));
            d+=dist(rng)/3;
            A.ge(j, j)=d;
            b.ge(j,0)=dist(rng)*dist(rng);
        }
        res=TomasSolver<dtype>::Solve(A, b);
        if (A*res!=b)
        {
            std::cerr<<"A"<<std::endl<<Matrix<dtype>(A)<<std::endl;
            std::cerr<<"b"<<std::endl<<b<<std::endl;
            std::cerr<<"b2"<<std::endl<<A*res<<std::endl;
            std::cerr<<"res"<<std::endl<<res<<std::endl;
            throw std::runtime_error("Something happend!");
        }
    }
}

std::size_t n,m;

int main(int argc, char* argv[])
{

    if (argc!=4)
    {
        std::cerr<<"wrong usage, need to have type, matrix size and number of test";
        return -1;
    }

    sscanf(argv[2], "%lu", &n);
    sscanf(argv[3], "%lu", &m);

    if (std::string(argv[1])=="d")
    {
        check<double>(n, m);
        return 0;
    }

    if (std::string(argv[1])=="ld") {
       check<long double>(n, m);
        return 0;
    }

    throw std::runtime_error("Bad parametrs");
}