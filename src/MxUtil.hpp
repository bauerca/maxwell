#ifndef MX_UTIL
#define MX_UTIL

#include <iostream>
#include <vector>
#include <cmath>

#include "MxTypes.h"
#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"

#include "Teuchos_XMLObject.hpp"

//#include "Epetra_CrsMatrix.h"

class Epetra_MultiVector;
class Epetra_CrsMatrix;
template<typename Scalar> class MxCrsMatrix;
class MxMap;

namespace MxUtil {

const double pi = 3.14159265358979323846;
const double lightspeed = 299792458.;
//const double dEps = std::numeric_limits<double>::epsilon();
const double dEps = 0.0;

template<typename T1, typename T2>
T1 toScalar(T2 val);

template<typename Scalar>
Scalar i();

template<>
inline
double i<double>() {return 0.0;}

template<>
inline
std::complex<double> i<std::complex<double> >() {
  return std::complex<double>(0.0, 1.0);
}

inline
int sign(double val) {
  if (val < -dEps) return -1;
  else if (val < dEps) return 0;
  else return 1;
}


template<typename T>
std::complex<T> convertScalar(T from, std::complex<T> & to) {
  to = std::complex<T>(from, 0);
  return to;
}

template<typename T>
T convertScalar(std::complex<T> from, T & to) {
  to = from.real();
  return to;
}

template<typename T>
T convertScalar(T from, T & to) {
  to = from;
  return to;
}

template<typename T1, typename T2>
void convertVector(std::vector<T1> const & v1, std::vector<T2> & v2) {
  size_t sz = v1.size();
  v2.resize(sz);
  for (size_t i = 0; i < sz; i++)
    v2[i] = T2(v1[i]);
}


#if 0
template<>
std::complex<double> toScalar<std::complex<double> >(double val) {
  return std::complex<double>(val, 0.0);
}

template<>
double toScalar<double>(double val) {
  return val;
}

template<>
std::complex<double> toScalar<std::complex<double> >(
std::complex<double> val) {
  return val;
}
#endif


int pow(int base, int exp);

namespace Strings {

template<typename T>
std::string typeToStr(T val) {
  std::stringstream out;
  out << val;
  return out.str();
}

template<typename T>
T strToType(const std::string & str);

template<>
int strToType<int>(const std::string & str);

template<>
std::string strToType<std::string>(const std::string & str);

template<>
double strToType<double>(const std::string & str);

template<typename T>
std::vector<T> strToStdVector(const std::string & str);

//template<>
//std::vector<std::string> strToStdVector<std::string>(const std::string & str);

std::string stripLeadingSpaces(std::string const & str);

std::string stripTrailingSpaces(std::string const & str);

std::string stripLTSpaces(std::string const & str);

std::string lhsEquation(std::string const & str);

std::pair<std::string, std::string> splitEquation(std::string const & str);

std::vector<std::string> lines(std::string const & str);

std::vector<std::string> split(std::string const & str, std::string sep);

}


namespace XML {

// this function gathers all lines within an xml object block that are not contained within
// a nested xml object block. Ignores blank lines.
std::vector<std::string> nodeLines(const ::Teuchos::XMLObject & node);

std::string getAttr(std::string name, const ::Teuchos::XMLObject & node, std::string dfault);

// no default return value. Returns an error if name is not found within XML block
std::string getAttr(std::string name, const ::Teuchos::XMLObject & node);

std::string getBodyAttr(std::string name, const ::Teuchos::XMLObject & node);

std::string getTagAttr(std::string name, const ::Teuchos::XMLObject & node);

}


template<typename T, size_t DIM>
MxDimMatrix<T, DIM> eye() {
  return MxDimMatrix<T, DIM>(MxDimVector<T, DIM>(ScalarTraits<T>::one()));
}



template<typename T>
inline
T mod(T x, T y) {return x - y * T(floor(double(x) / double(y)));}


template<typename T>
inline
void printStdVector(const std::vector<T> & in) {
  int size = in.size();
  std::cout << "[ ";
  for (int i = 0; i < size; i++) {
    if (i == size - 1)
      std::cout << in[i];
    else
      std::cout << in[i] << ", ";
  }
  std::cout << " ]\n";
}

template<typename T>
inline
bool inStdVector(T val, const std::vector<T> & v) {
  typename std::vector<T>::const_iterator iter;
  for (iter = v.begin(); iter != v.end(); ++iter)
    if (*iter == val) return true;
  return false;
}

// Type T should be in the family of unsigned ints
template<typename T>
std::vector<T> divisors(T val) {
  std::vector<T> res;
  for (int i = 2; i <= val; ++i) {
    if (val % i == 0) res.push_back(i);
  }
  return res;
}

template<typename T>
std::vector<std::vector<T> > unorderedFactorizations(T val, T maxFac) {
  std::vector<std::vector<T> > res;

  std::vector<T> div = divisors(val);
  //printStdVector(div);

  if (div.size() == 0) {
    res.push_back(div);
    return res;
  }
  else {
    std::vector<std::vector<T> > subfacs;

    typename std::vector<std::vector<T> >::const_iterator subfacIter;
    typename std::vector<T>::const_iterator iter1, iter2;
    for (iter1 = div.begin(); iter1 != div.end(); ++iter1) {
      if (*iter1 <= maxFac) {
        subfacs = unorderedFactorizations(val / (*iter1), *iter1);

        for (subfacIter = subfacs.begin(); subfacIter != subfacs.end(); ++subfacIter) {
          std::vector<T> fac(1, *iter1);
          for (iter2 = subfacIter->begin(); iter2 != subfacIter->end(); ++iter2) {
            fac.push_back(*iter2);
          }
          // remove the 1 from the back
          //if (fac[fac.size() - 1] == 1) fac.pop_back();
          res.push_back(fac);
          //printStdVector(fac);
        }
      }
    }
    return res;
  }
}

template<typename T>
std::vector<std::vector<T> > unorderedFactorizations(T val) {
  if (val == 1)
    return std::vector<std::vector<T> >(1, std::vector<T>(1, 1));
  else
    return unorderedFactorizations(val, val);
}



template<typename T>
inline
std::vector<std::vector<T> > permutations(std::vector<T> v) {

  std::vector<std::vector<T> > res, subsetPerms;
  std::vector<T> perm, subset, done;
  size_t size = v.size();
  //std::cout << size << "\n";
  if (size == 1)
    return std::vector<std::vector<T> >(1, v);
  else {
    T curVal;
    for (size_t i = 0; i < size; i++) {
      curVal = v[i];
      if (inStdVector(curVal, done)) continue;
      for (size_t j = i + 1; j < i + size; j++)
        subset.push_back(v[j % size]); 
      //printStdVector(subset);
      subsetPerms = permutations<T>(subset);
      subset.clear();

      for (size_t j = 0; j < subsetPerms.size(); j++) {
        perm.push_back(curVal);
        for (size_t k = 0; k < size - 1; k++)
          perm.push_back(subsetPerms[j][k]);
        res.push_back(perm);
        //printStdVector(perm);
        perm.clear();
      }
      done.push_back(curVal);
    }
    return res;
  }
}


template<typename T, int DIM>
MxDimVector<double, DIM> ndRootFind(const std::vector<const T *> & funcObjs, MXDIMVEC guess, double tol, int maxiter);


// templated over all objects that have methods 
// 'func', 'hasGradFunc', and possibly 'gradFunc'.
template<typename T, int DIM>
inline
MxDimVector<double, DIM> rootFind(const T & funcObj, MxDimVector<double, DIM> p1, MxDimVector<double, DIM> p2, double tol, int maxiter) {

  maxiter = 500;

  double len = (p2 - p1).norm();
  double scaleTol = tol * len;

  double f,df,f1,f2;
  double x1,xl,xh,tmp,dx,dxold,root;

  MxDimVector<double, DIM> dir = (p2 - p1) / len;
  double x2 = len;

  f1 = funcObj.func(p1);
  f2 = funcObj.func(p2);
  if (f1 == 0)
    return p1;
  if (f2 == 0)
    return p2;


  // put line search origin at p1
  x1=0.;
  if (f1<0) {
    xl = x1;
    xh = x2;
  }
  else {
    xl = x2;
    xh = x1;
  }
  root = 0.5*(xl + xh);
  dxold = fabs(xh - xl);
  dx = dxold;
  f = funcObj.func(p1 + dir * root);
  df = dir.dot(funcObj.gradFunc(p1 + dir * root)); // used to be negative
  for (int i=0; i<maxiter; i++) {
    if ((((root-xh)*df-f)*((root-xl)*df-f)>=0) || (fabs(2.*f)>fabs(dxold*df))) {
      dxold = dx;
      dx = 0.5*(xh-xl);
      root = xl + dx;
      if (xl==root) return p1 + dir*root;
    }
    else {
      dxold = dx;
      dx = f/df;
      tmp = root;
      root -= dx;
      if (tmp==root) return p1 + dir*root;
    }
    //std::cout << "root find acc: " << fabs(dx) << std::endl;
    if (fabs(dx) < scaleTol) return p1 + dir*root;
    f = funcObj.func(p1 + dir*root);
    df = dir.dot(funcObj.gradFunc(p1 + dir*root)); // used to be negative
    if (f<0) xl=root;
    else xh=root;
    //std::cout << root << std::endl << std::endl;
  }
  std::cout << "Root found to accuracy: " << fabs(dx) << std::endl;
  std::cout << "  p1 = "; p1.print();
  std::cout << "  f1 = " << funcObj.func(p1);
  std::cout << "  p2 = "; p2.print();
  std::cout << "  f2 = " << funcObj.func(p2);
  std::cout << "  dir = "; dir.print();
  return p1 + dir * root;
}

namespace Trilinos {

  template<typename Scalar>
  void massiveCrsMultiply(
    std::vector<RCP<MxCrsMatrix<Scalar> > > const & mats,
    std::vector<bool> const & transposes,
    //std::vector<char> transposes,
    const RCP<MxCrsMatrix<Scalar> > & res, bool fc = true);

  template<typename Scalar>
  RCP<MxCrsMatrix<Scalar> > identity(RCP<MxMap> map);
}

namespace Epetra {


  void stripZeros(Epetra_CrsMatrix & m);

  void removeConstField(Epetra_MultiVector & x);

  void printNorms(std::string name, const Epetra_MultiVector & x);

} // end namespace Epetra

namespace MxTpetra {

  template<typename Scalar>
  inline
  void getRowSums(Tpetra::CrsMatrix<Scalar, MxIndex> const & matrix,
      Tpetra::Vector<Scalar, MxIndex> & rowSums) {
    
    if (not rowSums.getMap()->isSameAs(*matrix.getRowMap())) {
      std::cout << "Tpetra::invRowSums : map mismatch.\n";
      throw 1;
    }

    Scalar sum = 0;
    size_t localInds = rowSums.getMap()->getNodeNumElements();
    Teuchos::ArrayView<Scalar> vals;
    Teuchos::ArrayView<MxIndex> indices;
    for (size_t i = 0; i < localInds; ++i) {
      matrix.getLocalRowView(i, indices, vals);

      sum = 0; 
      for (size_t j = 0; j < indices.size(); ++j)
        sum += vals[j];

      rowSums.replaceLocalValue(i, sum);

    }
  }

} // end namespace Tpetra

namespace HDF5 {

void saveArray(double *, int, const char *);

}

} // end namespace MxUtil

#endif
