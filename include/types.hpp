#ifndef types_hpp__
#define types_hpp__

#include <cstdint>
namespace LINARS {

template<typename dtype>
class IMatrix;

template<typename dtype>
class IIterator;

template<typename dtype>
class Matrix;

template<typename dtype>
class IMatrix;

template<typename dtype>
class Vector;

template<typename dtype>
class M3diag;

template<typename dtype>
class MCSR;

template<typename dtype>
class MDOK;

template<typename dtype>
class VMatrix;

template<typename dtype>
class MURtriang;

constexpr const uint32_t preset_max_iter=10000;
constexpr const long double preset_max_r=1e-9;

}

#endif