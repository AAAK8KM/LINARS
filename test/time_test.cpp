#include "imatrix.hpp"
#include "matrix.hpp"
#include "t2m.hpp"
#include "types.hpp"
#include <cstdint>
#include <fstream>
#include <ios>
#include <iterator>
#include <ostream>
#include <random>
#include <stdexcept>
#include <utility>
#include <iostream>
#include <vector>
#include <chrono>


std::random_device dev;
std::mt19937 rng(dev());
std::uniform_real_distribution<float> dist(-1000,1000);

const uint32_t N=3000;


std::vector<std::pair<uint32_t,uint32_t>> cs(N*N);

__attribute__((optimize("O0")))
int main()
{
    Vector<double> v(N), res(N),res2(N);
    for (uint32_t n=0;n<N*N;n++)
    {
        cs[n]=std::make_pair(n/N, n%N);
        v[n/N]=dist(rng);
    }
    std::fstream file("res.csv", std::ios_base::out);
    if (!file.is_open()) {
        std::cerr << "Error creating file"  << std::endl;
        return 1;
    }

    for (uint32_t n=0;n<=N*N;n+=N*N/20)
        for (uint32_t cnt=1;cnt<10;cnt++)
        {
            if (n==0) break;
            file<<n<<",";
            std::cout<<"prep "<<n<<" "<<cnt<<std::endl;
            std::shuffle(cs.begin(), cs.end(), rng);
            Matrix<double> A(std::make_pair(N, N));
            for (uint32_t i=0;i<n;i++)
                A.ge(cs[i].first, cs[i].second)=dist(rng);
            MCSR<double> A2(A);
            Matrix<double> A3(A);
            /*for (uint32_t i=0;i<n;i++)
                A3.ge(is[i], js[i])=dist(rng);
            MDOK<double> A4(A3);*/
            /*std::cout<<"ready"<<std::endl;
            std::cout<<A<<std::endl<<Matrix<double>(A2);
            std::cout<<Matrix<double>(res)<<std::endl<<Matrix<double>(A2*v);*/
            auto end = std::chrono::steady_clock::now();

            auto start = std::chrono::steady_clock::now();
            
            start = std::chrono::steady_clock::now();
            res=A2*v;
            end = std::chrono::steady_clock::now();
            
            if (res!=A3*v)
            {
                std::cout<<A<<std::endl<<Matrix<double>(A2);
                std::cout<<Matrix<double>(res)<<std::endl<<Matrix<double>(A2*v);
                
                throw std::runtime_error("Bad happens");
            }
            std::chrono::duration<double, std::milli> elapsed_ms = end - start;
            file<<elapsed_ms.count()<<",";
            
            start = std::chrono::steady_clock::now();
            res2=A*v;
            end = std::chrono::steady_clock::now();
            
            if (res2!=A3*v) throw std::runtime_error("Bad happens2");
            elapsed_ms = end - start;
            file<<elapsed_ms.count()<<std::endl;
        }
    file.flush();
    file.close();
}