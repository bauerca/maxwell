#ifndef MX_DIM_VEC
#define MX_DIM_VEC

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <string>
#include <complex>

// forward declaration because MxUtil uses MxDimVector.hpp
namespace MxUtil {
  namespace Strings {
template<typename T> std::vector<T> strToStdVector(const std::string & str);
  }
}


//#include "Matrix.h"
template<typename T, size_t DIM> class MxDimMatrix;

template<typename T, size_t DIM> 
class MxDimVector {
  public:
    MxDimVector();
    MxDimVector(T);
    
    template<size_t DIM2>
    explicit MxDimVector(const MxDimVector<T, DIM2> & v, T fill = 0) {
      for (size_t i = 0; i < DIM; ++i) {
        if (i < DIM2) data[i] = v[i];
        else data[i] = fill;
      }
    }

    T & operator[] (size_t i) {return data[i];}
    T operator[] (size_t i) const {return data[i];}

    MxDimVector<T, DIM> operator+ (const MxDimVector<T, DIM> &) const;
    MxDimVector<T, DIM> operator- (const MxDimVector<T, DIM> &) const;
    MxDimVector<T, DIM> operator* (const T &) const;
    MxDimVector<T, DIM> operator* (const MxDimVector<T, DIM> &) const;
    MxDimVector<T, DIM> operator* (const MxDimMatrix<T, DIM> &) const;
    MxDimVector<T, DIM> operator/ (const MxDimVector<T, DIM> &) const;
    MxDimVector<T, DIM> operator/ (const T &) const;
    MxDimVector<T, DIM> & operator/= (const T &);
    MxDimVector<T, DIM> & operator/= (const MxDimVector<T, DIM> &);
    MxDimVector<T, DIM> & operator-= (const MxDimVector<T, DIM> &);
    MxDimVector<T, DIM> & operator+= (const MxDimVector<T, DIM> &);
    MxDimVector<T, DIM> & operator*= (const T &);
    MxDimVector<T, DIM> & operator*= (const MxDimVector<T, DIM> &);
    bool operator== (const MxDimVector<T, DIM> &) const;
    bool operator!= (const MxDimVector<T, DIM> &) const;
    void print() const;
    T dot(const MxDimVector<T, DIM> &) const;
    T norm() const;
    T oneNorm() const;
    T sum() const;
    T prod() const;
    T max() const;
    T min() const;
    void fill(const T & val);
    void strFill(const std::string & str, T fill = 0);

  private:
    T data[DIM];

};

template<typename T, size_t DIM>
MxDimVector<T, DIM>::MxDimVector() {
  for (size_t i = 0; i < DIM; i++)
    data[i] = 0;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM>::MxDimVector(T val) {
  for (size_t i = 0; i < DIM; i++)
    data[i] = val;
}


template<typename T, size_t DIM>
MxDimVector<T, DIM> MxDimVector<T, DIM>::operator+(const MxDimVector<T, DIM> & v) const {
  MxDimVector<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    res[i] = data[i] + v[i];
  return res;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> operator-(const MxDimVector<T, DIM> & v) {
  MxDimVector<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    res[i] = -v[i];
  return res;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> MxDimVector<T, DIM>::operator-(const MxDimVector<T, DIM> & v) const {
  MxDimVector<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    res[i] = data[i] - v[i];
  return res;
}

template<typename T, size_t DIM>
T MxDimVector<T, DIM>::prod() const {
  T res(data[0]);
  for (size_t i = 1; i < DIM; i++)
    res *= data[i];
  return res;
}

template<typename T, size_t DIM>
T MxDimVector<T, DIM>::sum() const {
  T res(data[0]);
  for (size_t i = 1; i < DIM; i++)
    res += data[i];
  return res;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> MxDimVector<T, DIM>::operator* (const T & val) const {
  MxDimVector<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    res[i] = val * data[i];
  return res;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> MxDimVector<T, DIM>::operator/ (const T & val) const {
  MxDimVector<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    res[i] = data[i] / val;
  return res;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> MxDimVector<T, DIM>::operator/ (const MxDimVector<T, DIM> & v) const {
  MxDimVector<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    res[i] = data[i] / v[i];
  return res;
}


template<typename T, size_t DIM>
MxDimVector<T, DIM> MxDimVector<T, DIM>::operator* (const MxDimVector<T, DIM> & v) const {
  MxDimVector<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    res[i] = data[i] * v[i];
  return res;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> MxDimVector<T, DIM>::operator* (const MxDimMatrix<T, DIM> & m) const {
  MxDimVector<T, DIM> res(0);
  for (size_t i = 0; i < DIM; ++i)
    for (size_t j = 0; j < DIM; ++j)
      res[i] += data[j] * m(j, i);
  return res;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> & MxDimVector<T, DIM>::operator-= (const MxDimVector<T, DIM> & v) {
  for (size_t i = 0; i < DIM; i++)
    data[i] -= v[i];
  return *this;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> & MxDimVector<T, DIM>::operator+= (const MxDimVector<T, DIM> & v) {
  for (size_t i = 0; i < DIM; i++)
    data[i] += v[i];
  return *this;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> & MxDimVector<T, DIM>::operator/= (const MxDimVector<T, DIM> & v) {
  for (size_t i = 0; i < DIM; i++)
    data[i] /= v[i];
  return *this;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> & MxDimVector<T, DIM>::operator/= (const T & val) {
  for (size_t i = 0; i < DIM; i++)
    data[i] /= val;
  return *this;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> & MxDimVector<T, DIM>::operator*= (const T & val) {
  for (size_t i = 0; i < DIM; i++)
    data[i] *= val;
  return *this;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> & MxDimVector<T, DIM>::operator*= (const MxDimVector<T, DIM> & v) {
  for (size_t i = 0; i < DIM; i++)
    data[i] *= v[i];
  return *this;
}

template<typename T, size_t DIM>
bool MxDimVector<T, DIM>::operator==(const MxDimVector<T, DIM> & v) const {
  for (size_t i = 0; i < DIM; i++)
    if (data[i] != v[i]) return false;
  return true;
}

template<typename T, size_t DIM>
bool MxDimVector<T, DIM>::operator!=(const MxDimVector<T, DIM> & v) const {
  return !this->operator==(v);
}

template<typename T, size_t DIM>
T MxDimVector<T, DIM>::norm() const {
  return sqrt(this->dot(*this));
}

template<typename T, size_t DIM>
T MxDimVector<T, DIM>::oneNorm() const {
  T res = 0;
  for (size_t i = 0; i < DIM; ++i)
    res += fabs(data[i]);
  return res;
}

template<typename T, size_t DIM>
T MxDimVector<T, DIM>::dot(const MxDimVector<T, DIM> & v) const {
  T res = 0;
  for (size_t i = 0; i < DIM; i++)
    res += data[i] * v[i];
  return res;
}

template<typename T, size_t DIM>
void MxDimVector<T, DIM>::print() const {
  std::cout << "[ ";
  for (size_t i = 0; i < DIM; i++) std::cout << data[i] << ", ";
  std::cout<<" ]\n";
}

template<typename T, size_t DIM>
T MxDimVector<T, DIM>::max() const {
  T maxVal = data[0];
  for (size_t i = 0; i < DIM; i++)
    if (data[i] > maxVal)
      maxVal = data[i];
  return maxVal;
}

template<typename T, size_t DIM>
T MxDimVector<T, DIM>::min() const {
  T minVal = data[0];
  for (size_t i = 0; i < DIM; i++)
    if (data[i] < minVal)
      minVal = data[i];
  return minVal;
}

template<typename T, size_t DIM>
void MxDimVector<T, DIM>::fill(const T & val) {
  for (size_t i = 0; i < DIM; ++i)
    data[i] = val;
}

template<typename T, size_t DIM>
inline
void MxDimVector<T, DIM>::strFill(const std::string & str, T fill) {
  std::vector<T> strData = MxUtil::Strings::strToStdVector<T>(str);
  size_t size = strData.size(); 
  if (size < DIM)
    std::cout << "MxDimVector::strFill(...): Warning: not enough data in string, '"
              << str << "'. Padding with " << fill << ".\n";

  for (size_t i = 0; i < DIM; ++i) {
    if (i < size) data[i] = strData[i];
    else data[i] = fill;
  }
}





// non member function operators

template<typename T, size_t DIM>
MxDimVector<T, DIM> real(MxDimVector<std::complex<T>, DIM> v) {
  MxDimVector<T, DIM> res;
  for (size_t i = 0; i < DIM; ++i)
    res[i] = v[i].real();
  return res;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> imag(MxDimVector<std::complex<T>, DIM> v) {
  MxDimVector<T, DIM> res;
  for (size_t i = 0; i < DIM; ++i)
    res[i] = v[i].imag();
  return res;
}

template<typename T, size_t DIM>
std::ostream & operator<<(std::ostream & str, MxDimVector<T, DIM> const & v) {
  str << "[ ";
  for (size_t i = 0; i < DIM; i++) str << v[i] << ", ";
  str << "]\n";
  return str;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> operator* (const T & val, const MxDimVector<T, DIM> & v) {
  MxDimVector<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    res[i] = val * v[i];
  return res;
}

template<typename T, size_t DIM>
MxDimVector<T, DIM> operator/ (const T & val, const MxDimVector<T, DIM> & v) {
  MxDimVector<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    res[i] = val / v[i];
  return res;
}

template<typename T>
inline
MxDimVector<T, 3> cross(const MxDimVector<T, 2> & v1, const MxDimVector<T, 2> & v2) {
  MxDimVector<T, 3> res(0);
  res[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return res;
}

template<typename T>
inline
MxDimVector<T, 3> cross(const MxDimVector<T, 3> & v1, const MxDimVector<T, 3> & v2) {
  MxDimVector<T, 3> res;
  res[0] = v1[1] * v2[2] - v1[2] * v2[1];
  res[1] = v1[2] * v2[0] - v1[0] * v2[2];
  res[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return res;
}

template<typename T, size_t DIM>
inline
MxDimVector<T, DIM> fabs(const MxDimVector<T, DIM> & v) {
  MxDimVector<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    res[i] = fabs(v[i]);
  return res;
}

template<typename T, size_t DIM>
inline
bool anyEqual(MxDimVector<T, DIM> v1, MxDimVector<T, DIM> v2) {
  bool res = false;
  for (size_t i = 0; i < DIM; ++i)
    if (v1[i] == v2[i])
      res = true;
  return res;
}

// helper functions for easily creating 2,3-vectors
template<typename T>
MxDimVector<T, 2> vec2(T val1, T val2) {
  MxDimVector<T, 2> res;
  res[0] = val1; res[1] = val2;
  return res;
}

template<typename T>
MxDimVector<T, 3> vec3(T val1, T val2, T val3) {
  MxDimVector<T, 3> res;
  res[0] = val1; res[1] = val2; res[2] = val3;
  return res;
}

/*
// specialized constructors for 2,3-vectors
//

template<typename T> class MxDimVector<T, 2> {};
template<typename T> class MxDimVector<T, 3> {};



template<typename T>
MxDimVector<T, 2>::MxDimVector(T val1, T val2) {
  data[0] = val1; data[1] = val2;
}

template<typename T>
MxDimVector<T, 3>::MxDimVector(T val1, T val2, T val3) {
  data[0] = val1; data[1] = val2; data[2] = val3;
}
*/

#define MXDIMVEC MxDimVector<double, DIM>

typedef MxDimVector<double, 3> MxVecD3;
typedef MxDimVector<double, 2> MxVecD2;

#endif
