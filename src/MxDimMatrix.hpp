#ifndef MX_DIM_MATRIX
#define MX_DIM_MATRIX

#include <iostream>
#include <algorithm>
#include <complex>

#include "MxDimVector.hpp"
#include "Teuchos_LAPACK.hpp"
//#include "Teuchos_LAPACK_wrappers.hpp"

//extern "C" void dgetrf_(const int*, const int*, double*, const int*,
//                        const int*, int*);
//extern "C" void dgetri_(const int*, double*, const int*, const int*, 
//                        double*, const int*, int*);
//void dgesvd_(const char*, const char*, const int*, const int*, double*, const int*, double*, double*,
//             const int*, double*, const int*, double*, const int*, int*);



template<typename T, size_t DIM> 
class MxDimMatrix {
  public:
    // default ctor
    MxDimMatrix();

    // fill all elements with one value
    explicit MxDimMatrix(T);

    // make a diagonal matrix
    MxDimMatrix(const MxDimVector<T, DIM> &);

    // make an outer product matrix
    MxDimMatrix(const MxDimVector<T, DIM> &, const MxDimVector<T, DIM> &);

    static MxDimMatrix<T, DIM> I();

    MxDimVector<T, DIM> operator[] (size_t i) const;

    T & operator()(size_t i, size_t j) {return data[i][j];}
    T operator()(size_t i, size_t j) const {return data[i][j];}

    MxDimMatrix<T, DIM> operator+ (const MxDimMatrix<T, DIM> &) const;
    MxDimMatrix<T, DIM> operator- (const MxDimMatrix<T, DIM> &) const;
    MxDimMatrix<T, DIM> operator* (const T &) const;
    MxDimVector<T, DIM> operator* (const MxDimVector<T, DIM> &) const;
    MxDimMatrix<T, DIM> operator* (const MxDimMatrix<T, DIM> &) const;
    MxDimMatrix<T, DIM> operator/ (const T &) const;
    MxDimMatrix<T, DIM> & operator/= (const T &);
    MxDimMatrix<T, DIM> & operator-= (const MxDimMatrix<T, DIM> &);
    MxDimMatrix<T, DIM> & operator+= (const MxDimMatrix<T, DIM> &);
    MxDimMatrix<T, DIM> & operator*= (const T &);
    void setRow(size_t, const MxDimVector<T, DIM> &);
    void setCol(size_t, const MxDimVector<T, DIM> &);
    void getRow(size_t, MxDimVector<T, DIM> &) const;
    void getCol(size_t, MxDimVector<T, DIM> &) const;
    void print() const;
    double norm() const;
    T sum() const;
    T prod() const;
    T trace() const;
    //MxDimMatrix<T, DIM> & transpose();
    MxDimMatrix<T, DIM> transpose() const;
    MxDimMatrix<T, DIM> inv() const;
    T det() const;
    int eig(MxDimVector<std::complex<T>, DIM> & evals, MxDimMatrix<std::complex<T>, DIM> & levecs, MxDimMatrix<std::complex<T>, DIM> & revecs) const;
    int solve(MxDimVector<T, DIM> & x, const MxDimVector<T, DIM> & b) const;

  private:
    T data[DIM][DIM];

};

// default ctor: gotsta initialize yo
template<typename T, size_t DIM>
MxDimMatrix<T, DIM>::MxDimMatrix() {
  for (size_t i = 0; i < DIM; i++)
    for (size_t j = 0; j < DIM; j++)
      data[i][j] = 0;
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM>::MxDimMatrix(T val) {
  for (size_t i = 0; i < DIM; i++)
    for (size_t j = 0; j < DIM; j++)
      data[i][j] = val;
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM>::MxDimMatrix(const MxDimVector<T, DIM> & v) {
  for (size_t i = 0; i < DIM; i++) {
    for (size_t j = 0; j < DIM; j++) {
      if (i == j) data[i][j] = v[i];
      else data[i][j] = 0;
    }
  }
}

// outer product constructor
template<typename T, size_t DIM>
MxDimMatrix<T, DIM>::MxDimMatrix(const MxDimVector<T, DIM> & v1, const MxDimVector<T, DIM> & v2) {
  for (size_t i = 0; i < DIM; i++)
    for (size_t j = 0; j < DIM; j++)
      data[i][j] = v1[i] * v2[j];
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM> MxDimMatrix<T, DIM>::I() {
  MxDimMatrix<T, DIM> res(0);
  for (size_t i = 0; i < DIM; ++i)
    res(i, i) = 1;
  return res;
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM> MxDimMatrix<T, DIM>::operator+ (const MxDimMatrix<T, DIM> & m) const {
  MxDimMatrix<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    for (size_t j = 0; j < DIM; j++)
      res(i, j) = data[i][j] + m(i, j);
  return res;
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM> operator- (const MxDimMatrix<T, DIM> & m) {
  MxDimMatrix<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    for (size_t j = 0; j < DIM; j++)
      res(i, j) = -m(i, j);
  return res;
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM> MxDimMatrix<T, DIM>::operator- (const MxDimMatrix<T, DIM> & m) const {
  MxDimMatrix<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    for (size_t j = 0; j < DIM; j++)
      res(i, j) = data[i][j] - m(i, j);
  return res;
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM> MxDimMatrix<T, DIM>::operator* (const T & val) const {
  MxDimMatrix<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    for (size_t j = 0; j < DIM; j++)
      res(i, j) = val * data[i][j];
  return res;
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM> MxDimMatrix<T, DIM>::operator/ (const T & val) const {
  return this->operator*(T(1) / val);
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM> operator* (const T & val, const MxDimMatrix<T, DIM> & m) {
  MxDimMatrix<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    for (size_t j = 0; j < DIM; j++)
      res(i, j) = val * m(i, j);
  return res;
}

// matrix/vector multiplication
template<typename T, size_t DIM>
MxDimVector<T, DIM> MxDimMatrix<T, DIM>::operator* (const MxDimVector<T, DIM> & v) const {
  MxDimVector<T, DIM> res(0);
  for (size_t i = 0; i < DIM; i++)
    for (size_t j = 0; j < DIM; j++)
      res[i] += data[i][j] * v[j];
  return res;
}

// matrix/matrix multiplication
template<typename T, size_t DIM>
MxDimMatrix<T, DIM> MxDimMatrix<T, DIM>::operator* (const MxDimMatrix<T, DIM> & m) const {
  MxDimMatrix<T, DIM> res(0);
  for (size_t i = 0; i < DIM; i++)
    for (size_t j = 0; j < DIM; j++)
      for (size_t k = 0; k < DIM; k++)
        res(i, j) += data[i][k] * m(k, j);
  return res;
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM> & MxDimMatrix<T, DIM>::operator-= (const MxDimMatrix<T, DIM> & m) {
  for (size_t i = 0; i < DIM; i++)
    for (size_t j = 0; j < DIM; j++)
      data[i][j] -= m(i, j);
  return *this;
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM> & MxDimMatrix<T, DIM>::operator+= (const MxDimMatrix<T, DIM> & m) {
  for (size_t i = 0; i < DIM; i++)
    for (size_t j = 0; j < DIM; j++)
      data[i][j] += m(i, j);
  return *this;
}

//template<typename T, size_t DIM>
//MxDimMatrix<T, DIM> & MxDimMatrix<T, DIM>::transpose() {
//  for (size_t i = 0; i < DIM; i++)
//    for (size_t j = i + 1; j < DIM; j++)
//      std::swap(data[i][j], data[j][i]);
//  return *this;
//}
template<typename T, size_t DIM>
T MxDimMatrix<T, DIM>::trace() const {
  T res = 0;
  for (size_t i = 0; i < DIM; i++)
    res += data[i][i];
  return res;
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM> MxDimMatrix<T, DIM>::transpose() const {
  MxDimMatrix<T, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    for (size_t j = 0; j < DIM; j++)
      res(i, j) = data[j][i];
  return res;
}


template<typename T, size_t DIM>
MxDimVector<T, DIM> MxDimMatrix<T, DIM>::operator[] (size_t i) const {
  MxDimVector<T, DIM> res;
  for (size_t j = 0; j < DIM; ++j)
    res[j] = data[i][j];
  return res;
}

template<typename T, size_t DIM>
void MxDimMatrix<T, DIM>::print() const {
  for (size_t i = 0; i < DIM; i++)
    this->operator[](i).print();
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM> MxDimMatrix<T, DIM>::inv() const {

  MxDimMatrix<T, DIM> thisT = this->transpose(); 
  MxDimMatrix<T, DIM> solnT(MxDimVector<T, DIM>(1));

  int info;
  int ipiv[DIM];
  Teuchos::LAPACK<int, T> lapack;
  lapack.GESV(DIM, DIM, &thisT(0, 0), DIM, ipiv, &solnT(0, 0), DIM, &info);
  
  if (info != 0)
    std::cout << "MxDimMatrix::inv(): error in lapack routine. Return code: " << info << ".\n";

  return solnT.transpose();
}


template<typename T, size_t DIM>
int MxDimMatrix<T, DIM>::solve(MxDimVector<T, DIM> & x, const MxDimVector<T, DIM> & b) const {
  x = b;

  MxDimMatrix<T, DIM> copy(*this);

  Teuchos::LAPACK<int, T> lapack;
  MxDimVector<int, DIM> ipiv;
  //int ipiv[DIM];
  int info;

  lapack.GETRF(DIM, DIM, &copy(0, 0), DIM, &ipiv[0], &info);
  if (info != 0)
    std::cout << "MxDimMatrix::solve(...): error in lapack routine getrf. Return code: " << info << ".\n";

  lapack.GETRS('T', DIM, 1, &copy(0, 0), DIM, &ipiv[0], &x[0], DIM, &info);
  if (info != 0)
    std::cout << "MxDimMatrix::solve(...): error in lapack routine getrs. Return code: " << info << ".\n";

  return info;
}

template<typename T, size_t DIM>
int MxDimMatrix<T, DIM>::eig(MxDimVector<std::complex<T>, DIM> & evals,
MxDimMatrix<std::complex<T>, DIM> & levecs,
MxDimMatrix<std::complex<T>, DIM> & revecs) const {

  MxDimMatrix<T, DIM> thisT = this->transpose(); 

  MxDimVector<T, DIM> rva, iva;
  MxDimMatrix<T, DIM> rve, lve;

  int info;
  int lwork = 5 * DIM;
  T work[lwork];
  Teuchos::LAPACK<int, T> lapack;
  lapack.GEEV('V', 'V', DIM, &thisT(0, 0), DIM, &rva[0], &iva[0], &lve(0, 0), DIM, &rve(0, 0), DIM, work, lwork, &info);
  
  if (info != 0)
    std::cout << "MxDimMatrix::eig(): error in lapack routine. Return code: " << info << ".\n";

  // gather results
  for (size_t i = 0; i < DIM; ++i) {
    if (iva[i] != 0) {
      for (size_t k = 0; k < 2; ++k) {
        evals[i + k].real() = rva[i + k];
        evals[i + k].imag() = iva[i + k];
        for (size_t j = 0; j < DIM; ++j) {
          levecs(i + k, j).real() = lve(j, i);
          levecs(i + k, j).imag() = (k == 0 ? 1.0 : -1.0) * lve(j, i + 1);
          revecs(i + k, j).real() = rve(j, i);
          revecs(i + k, j).imag() = (k == 0 ? 1.0 : -1.0) * rve(j, i + 1);
        }
      }
      i++;
    }
    else {
      evals[i].real() = rva[i];
      evals[i].imag() = iva[i];
      for (size_t j = 0; j < DIM; ++j) {
        levecs(i, j).real() = lve(j, i);
        revecs(i, j).real() = rve(j, i);
      }
    }
  }

  return info;
}

template<typename T, size_t DIM>
inline
void MxDimMatrix<T, DIM>::setRow(size_t row, const MxDimVector<T, DIM> & v) {
  for (size_t i = 0; i < DIM; ++i)
    data[row][i] = v[i];
}

template<typename T, size_t DIM>
inline
void MxDimMatrix<T, DIM>::setCol(size_t col, const MxDimVector<T, DIM> & v) {
  for (size_t i = 0; i < DIM; ++i)
    data[i][col] = v[i];
}

template<typename T, size_t DIM>
inline
void MxDimMatrix<T, DIM>::getRow(size_t row, MxDimVector<T, DIM> & v) const {
  for (size_t i = 0; i < DIM; ++i)
    v[i] = data[row][i];
}

template<typename T, size_t DIM>
inline
void MxDimMatrix<T, DIM>::getCol(size_t col, MxDimVector<T, DIM> & v) const {
  for (size_t i = 0; i < DIM; ++i)
    v[i] = data[i][col];
}

// non member functions

template<typename T, size_t DIM>
std::ostream & operator<<(std::ostream & str, MxDimMatrix<T, DIM> const & m) {
  for (size_t i = 0; i < DIM; i++) str << m[i];
  return str;
}

// hard code the determinant since these are small matrices

template<typename T>
T det(MxDimMatrix<T, 1> const & m) {
  return m(0, 0);
}

template<typename T>
T det(MxDimMatrix<T, 2> const & m) {
  return m(0, 0)*m(1, 1) - m(0, 1)*m(1, 0);
}

template<typename T>
T det(MxDimMatrix<T, 3> const & m) {
  return m(0, 0) * m(1, 1) * m(2, 2) +
         m(0, 1) * m(1, 2) * m(2, 0) +
         m(0, 2) * m(1, 0) * m(2, 1) -
         m(0, 2) * m(1, 1) * m(2, 0) -
         m(0, 1) * m(1, 0) * m(2, 2) -
         m(0, 0) * m(1, 2) * m(2, 1);
}

// complex stuff
 
template<typename T, size_t DIM>
MxDimMatrix<T, DIM> real(MxDimMatrix<std::complex<T>, DIM> m) {
  MxDimMatrix<T, DIM> res;
  for (size_t i = 0; i < DIM; ++i)
    for (size_t j = 0; j < DIM; ++j)
      res(i, j) = m(i, j).real();
  return res;
}

template<typename T, size_t DIM>
MxDimMatrix<T, DIM> imag(MxDimMatrix<std::complex<T>, DIM> m) {
  MxDimMatrix<T, DIM> res;
  for (size_t i = 0; i < DIM; ++i)
    for (size_t j = 0; j < DIM; ++j)
      res(i, j) = m(i, j).imag();
  return res;
}


#endif
