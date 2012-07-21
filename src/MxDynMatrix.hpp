#ifndef MX_DYN_MATRIX
#define MX_DYN_MATRIX

#include <iostream>
#include <algorithm>
#include <vector>

#include "MxDimVector.hpp"

//#include "Teuchos_LAPACK_wrappers.hpp"

//extern "C" void dgetrf_(const int*, const int*, double*, const int*,
//                        const int*, int*);
//extern "C" void dgetri_(const int*, double*, const int*, const int*, 
//                        double*, const int*, int*);
//extern "C" void dgesvd_(const char*, const char*, const int*, const int*, 
//                        double*, const int*, double*, double*,
//                        const int*, double*, const int*, double*,
//                        const int*, int*);


template<typename T> 
class MxDynMatrix {
  public:
    MxDynMatrix() : rows(0), cols(0) {}

    // fill all elements with one value
    MxDynMatrix(size_t, size_t, const T = T());


    T & operator()(size_t i, size_t j) {return data[j + cols * i];}
    T operator()(size_t i, size_t j) const {return data[j + cols * i];}

    MxDynMatrix<T> operator+ (const MxDynMatrix<T> &) const;
    MxDynMatrix<T> operator- (const MxDynMatrix<T> &) const;
    MxDynMatrix<T> operator* (const T &) const;
    //MxDimVector<T> operator* (const MxDimVector<T> &) const;
    MxDynMatrix<T> operator* (const MxDynMatrix<T> &) const;
    MxDynMatrix<T> operator/ (const double &) const;
    MxDynMatrix<T> & operator/= (const double& dbl);
    MxDynMatrix<T> & operator-= (const MxDynMatrix<T> &);
    MxDynMatrix<T> & operator+= (const MxDynMatrix<T> &);
    MxDynMatrix<T> & operator*= (const T &);
    void print() const;
    double norm() const;
    T sum() const;
    T prod() const;
    //MxDynMatrix<T> & transpose();
    MxDynMatrix<T> transpose() const;
    size_t numRows() const {return rows;}
    size_t numCols() const {return cols;}


  private:
    size_t rows, cols;
    std::vector<T> data;

};

template<typename T>
MxDynMatrix<T>::MxDynMatrix(size_t nr, size_t nc, const T val) : 
rows(nr), cols(nc) {
  data.resize(rows * cols, val);
}

// matrix/matrix multiplication
template<typename T>
MxDynMatrix<T> MxDynMatrix<T>::operator* (const MxDynMatrix<T> & m) const {
  MxDynMatrix<T> res(rows, m.cols);
  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < m.cols; j++)
      for (size_t k = 0; k < cols; k++)
        res(i, j) += data[k + cols * i] * m.data[j + m.cols * k];
  return res;
}

//template<typename T>
//MxDynMatrix<T> & MxDynMatrix<T>::transpose() {
//  for (size_t i = 0; i < rows; i++)
//    for (size_t j = i + 1; j < cols; j++)
//      std::swap(data[j + cols * i], data[i + cols * j]);
//  return *this;
//}

template<typename T>
MxDynMatrix<T> MxDynMatrix<T>::transpose() const {
  MxDynMatrix<T> res(cols, rows);
  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      res(j, i) = data[j + i * cols];
  return res;
}

template<typename T>
void MxDynMatrix<T>::print() const {
  for (size_t i = 0; i < rows; i++) {
    std::cout << "[ ";
    for (size_t j = 0; j < cols; j++) {
      if (j == cols - 1)
        std::cout << this->operator()(i, j) << " ";
      else
        std::cout << this->operator()(i, j) << ", ";
    }
    std::cout << "]\n";
  }
}

//template<typename T>
//MxDynMatrix<T> operator* (const T & val, const MxDynMatrix<T> & m) {
//  MxDynMatrix<T> res;
//  for (size_t i = 0; i < DIM; i++)
//    for (size_t j = 0; j < DIM; j++)
//      res(i, j) = val * m(i, j);
//  return res;
//}
//
////template<typename T>
////MxDynMatrix<T>::MxDynMatrix(const MxDimVector<T> & v) {
////  for (size_t i = 0; i < DIM; i++) {
////    for (size_t j = 0; j < DIM; j++) {
////      if (i == j) data[i][j] = v[i];
////      else data[i][j] = 0;
////    }
////  }
////}
//
//template<typename T>
//MxDynMatrix<T> MxDynMatrix<T>::operator+ (const MxDynMatrix<T> & m) const {
//  MxDynMatrix<T> res;
//  for (size_t i = 0; i < DIM; i++)
//    for (size_t j = 0; j < DIM; j++)
//      res(i, j) = data[i][j] + m(i, j);
//  return res;
//}
//
//template<typename T>
//MxDynMatrix<T> operator- (const MxDynMatrix<T> & m) {
//  MxDynMatrix<T> res;
//  for (size_t i = 0; i < DIM; i++)
//    for (size_t j = 0; j < DIM; j++)
//      res(i, j) = -m(i, j);
//  return res;
//}
//
//template<typename T>
//MxDynMatrix<T> MxDynMatrix<T>::operator- (const MxDynMatrix<T> & m) const {
//  MxDynMatrix<T> res;
//  for (size_t i = 0; i < DIM; i++)
//    for (size_t j = 0; j < DIM; j++)
//      res(i, j) = data[i][j] - m(i, j);
//  return res;
//}
//
//template<typename T>
//MxDynMatrix<T> MxDynMatrix<T>::operator* (const T & val) const {
//  MxDynMatrix<T> res;
//  for (size_t i = 0; i < DIM; i++)
//    for (size_t j = 0; j < DIM; j++)
//      res(i, j) = val * data[i][j];
//  return res;
//}
//
//
//// matrix/vector multiplication
//template<typename T>
//MxDimVector<T> MxDynMatrix<T>::operator* (const MxDimVector<T> & v) const {
//  MxDimVector<T> res(0);
//  for (size_t i = 0; i < DIM; i++)
//    for (size_t j = 0; j < DIM; j++)
//      res[i] += data[i][j] * v[j];
//  return res;
//}
//
//// matrix/matrix multiplication
//template<typename T>
//MxDynMatrix<T> MxDynMatrix<T>::operator* (const MxDynMatrix<T> & m) const {
//  MxDynMatrix<T> res(0);
//  for (size_t i = 0; i < DIM; i++)
//    for (size_t j = 0; j < DIM; j++)
//      for (size_t k = 0; k < DIM; k++)
//        res(i, j) += data[i][k] * m(k, j);
//  return res;
//}
//
//template<typename T>
//MxDynMatrix<T> & MxDynMatrix<T>::operator-= (const MxDynMatrix<T> & m) {
//  for (size_t i = 0; i < DIM; i++)
//    for (size_t j = 0; j < DIM; j++)
//      data[i][j] -= m(i, j);
//  return *this;
//}
//
//template<typename T>
//MxDynMatrix<T> & MxDynMatrix<T>::operator+= (const MxDynMatrix<T> & m) {
//  for (size_t i = 0; i < DIM; i++)
//    for (size_t j = 0; j < DIM; j++)
//      data[i][j] += m(i, j);
//  return *this;
//}
//




#endif
