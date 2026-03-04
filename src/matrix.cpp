#define Dmatrix
#define Dvector

#include "imatrix.hpp"


template class LINARS::Matrix<double>;
template class LINARS::Matrix<float>;

template class LINARS::Vector<double>;
template class LINARS::Vector<float>;
double operator*(const LINARS::Vector<double>&,const LINARS::Vector<double>&);
float operator*(const LINARS::Vector<float>&,const LINARS::Vector<float>&);