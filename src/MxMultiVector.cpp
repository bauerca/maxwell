#include "MxMultiVector.hpp"

#include "MxUtil.hpp"
#include "MxMap.hpp"
#include "MxCrsMatrix.hpp"

#include "Epetra_Map.h"
#include "Tpetra_Map.hpp"


template<typename Scalar>
MxMultiVector<Scalar>::MxMultiVector(RCP<MxMap> map, size_t numVecs) : 
mMap(map),
mNumVecs(numVecs),
#if LINALG_BASE == TPETRA
mVec(Tpetra::createMultiVector<Scalar>(map->mMap, numVecs)) {}
#elif LINALG_BASE == EPETRA
mVec(rcp(new Epetra_MultiVector(
    ScalarTraits<Scalar>::isComplex ?
        *map->getComplexEpetraMap() : *map->mMap, int(numVecs)))) {
  if (ScalarTraits<Scalar>::isComplex)
    mEye = rcp(new MxCrsMatrix<Scalar>(
        mMap, MxUtil::i<Scalar>()));
}
#endif


template<typename Scalar>
MxMultiVector<Scalar>::MxMultiVector(MxMultiVector<Scalar> const & mv,
std::vector<size_t> const & vecInds, bool deepcopy) : 
mMap(mv.mMap),
mNumVecs(vecInds.size()),
#if LINALG_BASE == TPETRA
make error
#elif LINALG_BASE == EPETRA
mEye(mv.mEye) {
  int vInds[vecInds.size()];
  for (size_t i = 0; i < vecInds.size(); ++i) {
    vInds[i] = int(vecInds[i]);
  }
  mVec = rcp(new Epetra_MultiVector(deepcopy ? Copy : View,
    *mv.mVec, vInds, int(vecInds.size())));
}
#endif



// copy constructor
template<typename Scalar>
MxMultiVector<Scalar>::MxMultiVector(MxMultiVector<Scalar> const & mv) : 
mMap(mv.mMap),
mNumVecs(mv.getNumVecs()),
#if LINALG_BASE == TPETRA
mVec(rcp(new Tpetra::MultiVector<Scalar>(*mv.mVec))) {}
#elif LINALG_BASE == EPETRA
mVec(rcp(new Epetra_MultiVector(*mv.mVec))),
//mVec(rcp(new Epetra_MultiVector(Copy, *mv.mVec, 0, int(mv.getNumVecs())))),
mEye(mv.mEye) {}
#endif

// assignment operator
template<typename Scalar>
MxMultiVector<Scalar> & MxMultiVector<Scalar>::operator=(
MxMultiVector<Scalar> const & mv) {
  mMap = mv.mMap;
  mNumVecs = mv.mNumVecs;
#if LINALG_BASE == TPETRA
  mVec = rcp(new Tpetra::MultiVector<Scalar>(*mv.mVec));
#elif LINALG_BASE == EPETRA
  //mVec = rcp(new Epetra_MultiVector(*mv.mVec));

  (*mVec) = (*mv.mVec);

  //MxIndex len = this->getLocalLength();
  //for (size_t vec = 0; vec < mNumVecs;
  //for (MxIndex i = 0; i < len; ++i) {
  //  this->replaceLocalValue

  mEye = mv.mEye;
#endif
  return *this;
}


template<typename Scalar>
void MxMultiVector<Scalar>::set(Scalar val) {
  size_t n = mMap->getNodeNumIndices();
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < mNumVecs; ++j)
      replaceLocalValue(i, j, val);
}

// scale function

#if LINALG_BASE == TPETRA
template<typename Scalar>
void MxMultiVector<Scalar>::scale(Scalar val) {
  mVec->scale(val);
}
#elif LINALG_BASE == EPETRA
template<>
void MxMultiVector<double>::scale(double val) {
  mVec->Scale(val);
}

template<>
void MxMultiVector<MxComplex>::scale(MxComplex val) {
  size_t n = mMap->getNodeNumIndices();
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < mNumVecs; ++j)
      replaceLocalValue(i, j, val * this->operator()(i, j));
  
  //MxMultiVector<MxComplex> tmp(mMap, mNumVecs);
  //mEye->apply(*this, tmp);
  //mVec->Update(val.imag(), *tmp.mVec, val.real());
}
#endif


template<typename Scalar>
void MxMultiVector<Scalar>::scale(std::vector<Scalar> const & vals) {
  size_t n = mMap->getNodeNumIndices();
  for (size_t j = 0; j < mNumVecs; ++j)
    for (size_t i = 0; i < n; ++i)
      replaceLocalValue(i, j, vals[j] * (*this)(i, j));
}


#if LINALG_BASE == TPETRA
make error
#elif LINALG_BASE == EPETRA
template<>
void MxMultiVector<double>::conj() {}

template<>
void MxMultiVector<MxComplex>::conj() {
  size_t n = mMap->getNodeNumIndices();
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < mNumVecs; ++j)
      (*mVec)[j][2*i+1] *= -1.0;
}
#endif




template<typename Scalar>
void MxMultiVector<Scalar>::norm2(std::vector<double> & norms) const {
#if LINALG_BASE == TPETRA
  mVec->norm2(Teuchos::ArrayView<double>(norms));
#elif LINALG_BASE == EPETRA
  norms.resize(mNumVecs);
  mVec->Norm2(&norms[0]);
#endif
}

template<typename Scalar>
void MxMultiVector<Scalar>::normalize() {
  
  std::vector<double> norms;
  this->norm2(norms);

  std::vector<Scalar> invNorms(norms.size());
  Scalar dummy;
  for (size_t i = 0; i < mNumVecs; ++i)
    invNorms[i] = ScalarTraits<Scalar>::one() / 
      MxUtil::convertScalar(norms[i], dummy);
  
  size_t n = mMap->getNodeNumIndices();
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < mNumVecs; ++j)
      replaceLocalValue(i, j, this->operator()(i, j) / invNorms[j]);
}


#if LINALG_BASE == TPETRA
#elif LINALG_BASE == EPETRA
template<>
void MxMultiVector<double>::dot(MxMultiVector<double> const & mv,
std::vector<double> & res) const {
  res.resize(mNumVecs);
  mv.mVec->Dot(*mVec, &res[0]);
}

template<>
void MxMultiVector<MxComplex>::dot(MxMultiVector<MxComplex> const & mv,
std::vector<MxComplex> & res) const {
  res.resize(mNumVecs);

  // real part of dot products
  double resRe[mNumVecs];
  mVec->Dot(*mv.mVec, resRe);

  // imaginary part of dot products
  double resIm[mNumVecs];
  MxMultiVector<MxComplex> tmp(mMap, mNumVecs);
  mEye->apply(*this, tmp);
  tmp.mVec->Scale(-1.0);
  mVec->Dot(*tmp.mVec, resIm);

  for (size_t i = 0; i < mNumVecs; ++i)
    res[i] = MxComplex(resRe[i], resIm[i]);
}
#endif




#if LINALG_BASE == TPETRA
#elif LINALG_BASE == EPETRA
template<>
void MxMultiVector<double>::update(double scalarA, MxMultiVector<double> const & mvA,
double scalarThis) {
  mVec->Update(scalarA, *mvA.mVec, scalarThis);
}

template<>
void MxMultiVector<MxComplex>::update(MxComplex scalarA,
MxMultiVector<MxComplex> const & mvA, MxComplex scalarThis) {

  // copy to tmp first in case mvA == *this
  MxMultiVector<MxComplex> tmp(mvA);
  this->scale(scalarThis);
  tmp.scale(scalarA);

  mVec->Update(1.0, *tmp.mVec, 1.0);
}
#endif

#if 0
#if LINALG_BASE == TPETRA
#elif LINALG_BASE == EPETRA
template<>
void MxMultiVector<double>::multiply(
    char transA, char transB, double scalarAB,
    MxMultiVector<double> const & mvA,
    MxMultiVector<double> const & mvB,
    double scalarThis) {
  mVec->Multiply(transA, transB, scalarAB,
      *mvA.mVec, *mvB.mVec, scalarThis);
}

template<>
void MxMultiVector<MxComplex>::multiply(
    char transA, char transB, MxComplex scalarAB,
    MxMultiVector<MxComplex> const & mvA,
    MxMultiVector<MxComplex> const & mvB,
    MxComplex scalarThis) {

  this->scale(scalarThis);
  MxMultiVector<MxComplex> tmp(mvA);
  tmp.scale(scalarA);

  mVec->Update(1.0, *tmp.mVec, 1.0);
}
#endif
#endif


#if LINALG_BASE == TPETRA
#elif LINALG_BASE == EPETRA
template<typename Scalar>
void MxMultiVector<Scalar>::setVector(size_t vecIndex,
RCP<MxVector<Scalar> const> vec, bool deepcopy) {
  if (deepcopy) {
    for (int i = 0; i < mVec->Map().NumMyElements(); ++i)
      mVec->ReplaceMyValue(i, int(vecIndex), (*vec->mVec)[0][i]);
  }
  else
    (*mVec)(vecIndex) = (*vec->mVec)(0);
}
#endif

#if LINALG_BASE == TPETRA
#elif LINALG_BASE == EPETRA
template<typename Scalar>
RCP<MxVector<Scalar> const> MxMultiVector<Scalar>::getVector(
size_t vecIndex, bool deepcopy) const {
  RCP<MxVector<Scalar> > res =
      rcp(new MxVector<Scalar>(*this, vecIndex, deepcopy));
  return res;
}
#endif

#if LINALG_BASE == TPETRA
#elif LINALG_BASE == EPETRA
template<typename Scalar>
RCP<MxVector<Scalar> > MxMultiVector<Scalar>::getVectorNonConst(
size_t vecIndex, bool deepcopy) {
  RCP<MxVector<Scalar> > res =
      rcp(new MxVector<Scalar>(*this, vecIndex, deepcopy));
  return res;
}
#endif

#if 0
template<typename Scalar>
RCP<Epetra_MultiVector> MxMultiVector<Scalar>::getEpetraMultiVector() {
  if (mEpetraVec == Teuchos::null) {
    RCP<Epetra_Map> map(mMap->getEpetraMap<Scalar>());

    mEpetraVec = rcp(new Epetra_MultiVector(*map, int(mNumVecs)));

    size_t numInds = mMap->getNodeNumIndices();
    bool isComplex = ScalarTraits<Scalar>::isComplex;

    Teuchos::ArrayRCP<const Scalar> vals;
    for (size_t i = 0; i < mNumVecs; ++i) {
      vals = mVec->getData(i);
      for (size_t j = 0; j < numInds; ++j) {
        if (isComplex) {
          mEpetraVec->ReplaceMyValue(2*j, i,
              ScalarTraits<Scalar>::real(vals[j]));
          mEpetraVec->ReplaceMyValue(2*j + 1, i,
              ScalarTraits<Scalar>::imag(vals[j]));
        }
        else
          mEpetraVec->ReplaceMyValue(j, i,
              ScalarTraits<Scalar>::real(vals[j]));
      }
    }
  }
  return mEpetraVec;
}
#endif

template class MxMultiVector<double>;
template class MxMultiVector<MxComplex>;
