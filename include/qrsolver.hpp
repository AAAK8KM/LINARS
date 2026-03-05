#ifndef qrsolver_h__
#define qrsolver_h__

#include "imatrix.hpp"
#include "qrdec.hpp"
#include "gauss.hpp"
#include "t2m.hpp"

namespace LINARS{

template <typename dtype,typename  Mtype>
requires IsMatrix<dtype, Mtype>
using QRDretT = std::invoke_result_t<decltype(QRdecompositionH<dtype,Mtype>), const Mtype&>(const Mtype&);

template <typename dtype,typename  Mtype>
requires IsMatrix<dtype, Mtype>
using QRSretT = std::conditional_t<std::is_same<Mtype, Vector<dtype>>::value, Vector<dtype>, VMatrix<dtype>>;

template <typename dtype,typename  Mtype1,typename  Mtype2>
requires IsMatrix<dtype, Mtype1> && IsMatrix<dtype, Mtype2>
QRSretT<dtype,Mtype2> QRSolver(const Mtype1& A, const Mtype2& b, QRDretT<dtype,Mtype1>& QRdec_f)
{
    auto [Q,R] = QRdec_f(A);
    Q.reverse();
    return ReverseGauss(R,Q(b));
}

}

#endif