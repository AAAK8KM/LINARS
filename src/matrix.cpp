#define Dmatrix
#define Dvector

#include "imatrix.hpp"

namespace LINARS{

template class Matrix<double>;
template class Matrix<float>;
template class VMatrix<double>;
template class VMatrix<float>;

template class Vector<double>;
template class Vector<float>;
//double operator*(const Vector<double>&,const Vector<double>&);
//float operator*(const Vector<float>&,const Vector<float>&);

}