#include "fstream"
#include "matrix.hpp"
#include "t2m.hpp"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>

std::size_t n,m;

template <typename T>
void check(std::fstream& f)
{
        f>>n>>m;
        Matrix<T> A(n,m);
        f>>A;
        f>>n>>m;
        Matrix<T> B(n,m);
        f>>B;
        f>>n>>m;
        Matrix<T> C(n,m);
        f>>C;
        if (A*B!=C)
            throw std::runtime_error("FAIL");
}

int main(int argc, char* argv[])
{

    if (argc!=3)
    {
        std::cerr<<"wrong usage, need to have type and path to sample";
        return -1;
    }

    std::fstream f;
    f.open((std::string(argv[2])));

    if (std::string(argv[1])=="d")
    {
        check<double>(f);
        return 0;
    }

    if (std::string(argv[1])=="f")
    {
        check<float>(f);
        return 0;
    }
}