#include "fstream"
#include "matrix.hpp"
#include "t2m.hpp"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>

std::size_t n,m;

using namespace LINARS;

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
        double d;
        f>>n>>m;
        Matrix<double> M(n,m);
        f>>M;
        f>>n>>m;
        for (std::size_t i=0;i<n;i++)
            for (std::size_t j=0;j<m;j++)
                {
                    f>>d;
                    if (d!=M.ge(i, j))
                        throw std::runtime_error("FAIL");
                }
        return 0;
    }

    if (std::string(argv[1])=="f")
    {
        float d;
        f>>n>>m;
        Matrix<float> M(n,m);
        f>>M;
        f>>n>>m;
        for (std::size_t i=0;i<n;i++)
            for (std::size_t j=0;j<m;j++)
                {
                    f>>d;
                    if (d!=M.ge(i, j))
                        throw std::runtime_error("FAIL");
                }
        return 0;
    }
}