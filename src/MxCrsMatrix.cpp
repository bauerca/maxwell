#include "MxCrsMatrix.hpp"

#include "MxMap.hpp"
#include "MxUtil.hpp"

#include "Epetra_CrsMatrix.h"

#include "TpetraExt_MatrixMatrix.hpp"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

template<typename Scalar>
MxCrsMatrix<Scalar>::MxCrsMatrix(RCP<MxMap> rowMap) : 
mFilled(false),
mRowMap(rowMap),
mColMap(Teuchos::null) {
  initBaseMatrix();
}

template<typename Scalar>
MxCrsMatrix<Scalar>::MxCrsMatrix(RCP<MxMap> rowMap, RCP<MxMap> colMap) : 
mFilled(false),
mRowMap(rowMap),
mColMap(colMap) {
  initBaseMatrix();
}

template<typename Scalar>
MxCrsMatrix<Scalar>::MxCrsMatrix(RCP<MxMap> rowMap, Scalar diag) : 
mFilled(false),
mRowMap(rowMap),
mColMap(rowMap) {
  initBaseMatrix();

  MxIndex row;
  size_t n = rowMap->getNodeNumIndices();
  for (size_t i = 0; i < n; ++i) {
    row = rowMap->getGlobalIndex(i);
    insertRowValues(row, 1, &row, &diag);
  }
  fillComplete(rowMap, rowMap);
}

template<typename Scalar>
MxCrsMatrix<Scalar>::MxCrsMatrix(RCP<MxVector<Scalar> > diag) : 
mFilled(false),
mRowMap(diag->mMap),
mColMap(diag->mMap) {
  initBaseMatrix();

  MxIndex row;
  Scalar val;
  size_t n = mRowMap->getNodeNumIndices();
  for (size_t i = 0; i < n; ++i) {
    row = mRowMap->getGlobalIndex(i);
    val = diag->getValue(i);
    insertRowValues(row, 1, &row, &val);
  }
  fillComplete(mRowMap, mRowMap);
}

template<typename Scalar>
void MxCrsMatrix<Scalar>::initBaseMatrix() {
#if LINALG_BASE == TPETRA
  mMatrix = rcp(new Tpetra::CrsMatrix<Scalar, MxIndex>(mRowMap->mMap, 0));
#elif LINALG_BASE == EPETRA
  if (mColMap == Teuchos::null) {
    mMatrix = rcp(new Epetra_CrsMatrix(Copy,
      ScalarTraits<Scalar>::isComplex ?
        *mRowMap->getComplexEpetraMap() : *mRowMap->mMap, 0));
  } else {
    mMatrix = rcp(new Epetra_CrsMatrix(Copy,
      ScalarTraits<Scalar>::isComplex ?
        *mRowMap->getComplexEpetraMap() : *mRowMap->mMap,
      ScalarTraits<Scalar>::isComplex ?
        *mColMap->getComplexEpetraMap() : *mColMap->mMap,
      0));
  }
  //std::cout << mRowMap->getComplexEpetraMap()->NumGlobalElements() << "\n";
#endif
}


#if LINALG_BASE == TPETRA
template<typename Scalar>
void MxCrsMatrix<Scalar>::insertRowValues(MxIndex row, size_t numEntries,
MxIndex const * cols, Scalar const * vals) {
  mMatrix->insertGlobalValues(row,
    Teuchos::ArrayView<MxIndex const>(cols, numEntries),
    Teuchos::ArrayView<Scalar const>(vals, numEntries));
}
#elif LINALG_BASE == EPETRA
template<>
void MxCrsMatrix<double>::insertRowValues(MxIndex row,
size_t numEntries, MxIndex const * cols, double const * vals) {
  int colsInt[numEntries];
  for (size_t i = 0; i < numEntries; ++i)
    colsInt[i] = int(cols[i]);
  int res = mMatrix->InsertGlobalValues(int(row), int(numEntries), vals,
    colsInt);

  if (res < 0) {
    std::cout << "MxCrsMatrix::insertRowValues: Epetra error code " << res << ". row=" << row <<
      " numEntries=" << numEntries << "\n";
    throw 1;
  }
}

template<>
void MxCrsMatrix<MxComplex>::insertRowValues(MxIndex row,
size_t numEntries, MxIndex const * inds, MxComplex const * vals) {
  size_t n = numEntries;

  std::vector<double> valsRe(n);
  std::vector<double> valsIm(n);
  std::vector<double> valsImNeg(n);

  std::vector<int> indsRe(n);
  std::vector<int> indsIm(n);

  for (size_t i = 0; i < n; ++i) {
    valsRe[i] = vals[i].real();
    valsIm[i] = vals[i].imag();
    valsImNeg[i] = -vals[i].imag();
    indsRe[i] = 2*inds[i];
    indsIm[i] = 2*inds[i] + 1;
  }

  mMatrix->InsertGlobalValues(2*row, n, &valsRe[0], &indsRe[0]);
  mMatrix->InsertGlobalValues(2*row, n, &valsImNeg[0], &indsIm[0]);
  mMatrix->InsertGlobalValues(2*row+1, n, &valsIm[0], &indsRe[0]);
  mMatrix->InsertGlobalValues(2*row+1, n, &valsRe[0], &indsIm[0]);

}
#endif


template<typename Scalar>
void MxCrsMatrix<Scalar>::insertRowValues(MxIndex row,
std::vector<MxIndex> cols, std::vector<Scalar> vals) {
  insertRowValues(row, cols.size(), &cols[0], &vals[0]);
}


#if LINALG_BASE == TPETRA
template<typename Scalar>
void MxCrsMatrix<Scalar>::getLocalRow(MxIndex localIndex, std::vector<MxIndex> & localIndices,
std::vector<Scalar> & vals) const {
  size_t rowEntries = mMatrix->getNumEntriesInLocalRow(localIndex);
  localIndices.resize(rowEntries);
  vals.resize(rowEntries);

  mMatrix->getLocalRowCopy(localIndex,
    Teuchos::ArrayView<MxIndex>(localIndices),
    Teuchos::ArrayView<Scalar>(vals));
}
#elif LINALG_BASE == EPETRA
template<>
void MxCrsMatrix<double>::getLocalRow(MxIndex localIndex, std::vector<MxIndex> & localIndices,
std::vector<double> & vals) const {
  int rowEntries = mMatrix->NumMyEntries(localIndex);
  localIndices.resize(rowEntries);
  vals.resize(rowEntries);

  int dummy;
  int intInds[rowEntries];
  mMatrix->ExtractMyRowCopy(localIndex, rowEntries, dummy, &vals[0], intInds);

  for (size_t i = 0; i < rowEntries; ++i)
    localIndices[i] = MxIndex(intInds[i]);
}

template<>
void MxCrsMatrix<MxComplex>::getLocalRow(MxIndex localIndex, std::vector<MxIndex> & localIndices,
std::vector<MxComplex> & vals) const {
  int rowEntries = mMatrix->NumMyEntries(2*localIndex);
  std::vector<int> inds(rowEntries);
  std::vector<double> vls(rowEntries);
  if (rowEntries % 2 != 0)
    std::cout << "Problem in MxCrsMatrix<Complex>::getLocalRow: rowEntries % 2 != 0\n";
  localIndices.resize(rowEntries / 2);
  vals.resize(rowEntries / 2);

  int dummy;
  mMatrix->ExtractMyRowCopy(2*localIndex, rowEntries, dummy, &vls[0], &inds[0]);

  int indRe, indIm;
  size_t j = 0;
  for (int i = 0; i < rowEntries; i += 2) {
    indRe = inds[i];
    indIm = inds[i+1];

    if (indRe % 2 != 0 or indIm % 2 != 1)
      std::cout << "Problem in MxCrsMatrix<Complex>::getLocalRow: extracted indices\n";

    vals[j] = MxComplex(vls[i], -vls[i+1]);
    localIndices[j] = indRe / 2;
    j++;
  }
}
#endif



#if LINALG_BASE == TPETRA
template<typename Scalar>
void MxCrsMatrix<Scalar>::scale(Scalar val) {
  mMatrix->scale(val);
}
#elif LINALG_BASE == EPETRA
template<>
void MxCrsMatrix<double>::scale(double val) {
  mMatrix->Scale(val);
}

template<>
void MxCrsMatrix<MxComplex>::scale(MxComplex val) {
  // make diagonal matrix
  MxCrsMatrix<MxComplex> valM(mRowMap, val);

  // resultant matrix
  MxCrsMatrix<MxComplex> valXthis(mRowMap, mColMap);

  // res = scale * this
  MxCrsMatrix<MxComplex>::matrixMatrixMultiply(valM, false, *this, false,
    valXthis, mFilled);

  // replace
  mMatrix = valXthis.mMatrix;
}
#endif






#if LINALG_BASE == TPETRA
template<typename Scalar>
void MxCrsMatrix<Scalar>::leftScale(RCP<MxVector<Scalar> > vec) {
  mMatrix->leftScale(*vec->mVec->getVectorNonConst(0));
}
#elif LINALG_BASE == EPETRA
template<typename Scalar>
void MxCrsMatrix<Scalar>::leftScale(RCP<MxVector<Scalar> > vec) {
  MxCrsMatrix<Scalar> vecOp(vec);
  MxCrsMatrix<Scalar> res(mRowMap, mColMap);

  // res = scale * this
  MxCrsMatrix<Scalar>::matrixMatrixMultiply(vecOp, false, *this, false,
    res, mFilled);

  // reset underlying matrix of this to newly created underlying matrix
  mMatrix = res.mMatrix;
  // res object deleted, and old this->mMatrix deleted
}
#endif



// generic to any linalg_base
template<typename Scalar>
RCP<MxVector<Scalar> > MxCrsMatrix<Scalar>::getRowSums(bool inverse) const {
  RCP<MxVector<Scalar> > res;
  res = rcp(new MxVector<Scalar>(mRowMap));

  size_t n = mRowMap->getNodeNumIndices();
  std::vector<Scalar> vals;
  std::vector<MxIndex> inds;
  Scalar sum, dummy;
  for (size_t i = 0; i < n; ++i) {
    getLocalRow(i, inds, vals);
    for (size_t j = 0; j < inds.size(); ++j)
      sum += MxUtil::convertScalar(ScalarTraits<Scalar>::magnitude(vals[j]), dummy);

    if (inverse)
      res->replaceLocalValue(i, ScalarTraits<Scalar>::one() / sum);
    else
      res->replaceLocalValue(i, sum);
  }
 
  return res;
}




template<typename Scalar>
void MxCrsMatrix<Scalar>::fillComplete(RCP<MxMap> domainMap,
RCP<MxMap> rangeMap) {
#if LINALG_BASE == TPETRA
  mMatrix->fillComplete(domainMap->mMap, rangeMap->mMap);
  //mColMap = rcp(new MxMap(mMatrix->getColMap(), mRowMap->getComm()));
#elif LINALG_BASE == EPETRA
  //std::cout << "Starting fill complete.\n";
  mColMap = domainMap;
  mMatrix->FillComplete(
      ScalarTraits<Scalar>::isComplex ?
          *domainMap->getComplexEpetraMap() : *domainMap->mMap,
      ScalarTraits<Scalar>::isComplex ?
          *rangeMap->getComplexEpetraMap() : *rangeMap->mMap);
  //std::cout << "Fill complete finished.\n";
  //mColMap = rcp(new MxMap(mMatrix->getColMap(), mRowMap->getComm()));
#endif
  mFilled = true;
}



template<typename Scalar>
void MxCrsMatrix<Scalar>::apply(MxMultiVector<Scalar> const & x,
MxMultiVector<Scalar> & y) const {
#if LINALG_BASE == TPETRA
#elif LINALG_BASE == EPETRA
  mMatrix->Apply(*x.mVec, *y.mVec);
#endif
}






template<typename Scalar>
void MxCrsMatrix<Scalar>::matrixMatrixMultiply(
MxCrsMatrix<Scalar> const & a, bool transA,
MxCrsMatrix<Scalar> const & b, bool transB,
MxCrsMatrix<Scalar> & result, bool fc) {
#if LINALG_BASE == TPETRA
  Tpetra::MatrixMatrix::Multiply(*a.mMatrix, transA,
    *b.mMatrix, transB, *result.mMatrix, fc);
#elif LINALG_BASE == EPETRA
  EpetraExt::MatrixMatrix::Multiply(*a.mMatrix, transA,
    *b.mMatrix, transB, *result.mMatrix, fc);
#endif
  if (fc) {
    result.mFilled = true;
    if (transB)
      result.mColMap = b.mRowMap;
    else
      result.mColMap = b.mColMap;

    // don't set result row map because user should have
    // already set that.
  }
}





#if LINALG_BASE == TPETRA
template<typename Scalar>
void MxCrsMatrix<Scalar>::matrixMatrixAdd(
MxCrsMatrix<Scalar> const & a, bool transA, Scalar scalarA,
MxCrsMatrix<Scalar> const & b, bool transB, Scalar scalarB,
MxCrsMatrix<Scalar> & result, bool fc) {
  Tpetra::MatrixMatrix::Add(*a.mMatrix, transA, scalarA,
    *b.mMatrix, transB, scalarB, result.mMatrix);
  if (fc) {
    if (transA)
      result.fillComplete(a.getDomainMap(), a.getRangeMap());
    else
      result.fillComplete(a.getRangeMap(), a.getDomainMap());
  }
}
#elif LINALG_BASE == EPETRA
template<>
void MxCrsMatrix<double>::matrixMatrixAdd(
MxCrsMatrix<double> const & a, bool transA, double scalarA,
MxCrsMatrix<double> const & b, bool transB, double scalarB,
MxCrsMatrix<double> & result, bool fc) {
  Epetra_CrsMatrix * cPtr = result.mMatrix.get();
  EpetraExt::MatrixMatrix::Add(*a.mMatrix, transA, scalarA,
    *b.mMatrix, transB, scalarB, cPtr);
  if (fc) {
    if (a.isFilled()) {
      if (transA)
        result.fillComplete(a.getDomainMap(), a.getRangeMap());
      else
        result.fillComplete(a.getRangeMap(), a.getDomainMap());
    } else if (b.isFilled()) {
      if (transB)
        result.fillComplete(b.getDomainMap(), b.getRangeMap());
      else
        result.fillComplete(b.getRangeMap(), b.getDomainMap());
    } else {
      std::cout << "a nor b is filled\n";
      throw 1;
    }
  }
}

template<>
void MxCrsMatrix<MxComplex>::matrixMatrixAdd(
MxCrsMatrix<MxComplex> const & a, bool transA, MxComplex scalarA,
MxCrsMatrix<MxComplex> const & b, bool transB, MxComplex scalarB,
MxCrsMatrix<MxComplex> & result, bool fc) {
  
  // make diagonal matrices
  MxCrsMatrix<MxComplex> sA(a.getRangeMap(), scalarA);
  MxCrsMatrix<MxComplex> sB(b.getRangeMap(), scalarB);

  MxCrsMatrix<MxComplex> sAxA(a.getRangeMap(), a.getDomainMap());
  MxCrsMatrix<MxComplex> sBxB(b.getRangeMap(), b.getDomainMap());

  MxCrsMatrix<MxComplex>::matrixMatrixMultiply(sA, false, a, false, sAxA, true);
  MxCrsMatrix<MxComplex>::matrixMatrixMultiply(sB, false, b, false, sBxB, true);

  // does not call fill complete
  Epetra_CrsMatrix * cPtr = result.mMatrix.get();
  EpetraExt::MatrixMatrix::Add(*sAxA.mMatrix, transA, 1.0,
    *sBxB.mMatrix, transB, 1.0, cPtr);
  if (fc) {
    if (transA)
      result.fillComplete(a.getDomainMap(), a.getRangeMap());
    else
      result.fillComplete(a.getRangeMap(), a.getDomainMap());
  }
}
#endif


template<typename Scalar>
void MxCrsMatrix<Scalar>::save(std::string name) const {
#if LINALG_BASE == EPETRA
  std::string nname = name + std::string(".mm");
  EpetraExt::RowMatrixToMatrixMarketFile(nname.c_str(), *mMatrix);
#endif
}


#if 0
template<>
RCP<Tpetra::CrsMatrix<double, MxIndex> >
MxCrsMatrix<double>::getTpetraCrsMatrix() const {
  RCP<Tpetra::CrsMatrix<double, MxIndex> > res;
#if CRS_BASE == TPETRA
  res = mMatrix;
#elif CRS_BASE == EPETRA
  size_t numRowElems = mRowMap.NumMyElements();
  MxIndex * rowElems = mRowMap.MyGlobalElements();
  MxIndex * colElems = mColMap.MyGlobalElements();

  std::vector<int> globColInds(mMatrix->MaxNumEntries());
  
  Teuchos::ArrayView<MxIndex> rowElemsView(rowElems, numRowElems);
  res = Tpetra::createCrsMatrix(Tpetra::createNonContigMap(rowElemsView));

  // now fill the matrix
  int numEntries;
  int * localColInds;
  double * vals;
  for (size_t i = 0; i < numRowElems; ++i) {
    mMatrix->ExtractMyRowView(i, numEntries, vals, localColInds);

    for (size_t j = 0; j < numEntries; ++j)
      globColInds[j] = colElems[localColInds[j]];

    res->insertGlobalValues(rowElems[i],
      Teuchos::ArrayView<MxIndex>(&globColInds[0], numEntries);
      Teuchos::ArrayView<double>(vals, numEntries));
  }
  res->fillComplete(mColMap, mRowMap);
#endif
  return res;
}



template<typename Scalar>
RCP<Epetra_CrsMatrix> MxCrsMatrix<Scalar>::getEpetraCrsMatrix() const {

  RCP<Epetra_Map> rowMap = mRowMap->getEpetraMap<Scalar>();
  RCP<Epetra_Map> colMap = mColMap->getEpetraMap<Scalar>();

  RCP<Epetra_CrsMatrix> res = rcp(new Epetra_CrsMatrix(Copy,
    *rowMap, 0));

  size_t numRowInds = mRowMap->getNodeNumIndices();
  MxIndex const * rowInds = mRowMap->getNodeIndexList();
  MxIndex const * colInds = mColMap->getNodeIndexList();

  std::vector<MxIndex> globColInds(mMatrix->getNodeMaxNumRowEntries());

  // now fill the matrix
  int numEntries;
  Teuchos::ArrayView<MxIndex const> localColInds;
  Teuchos::ArrayView<Scalar const> vals;
  for (size_t i = 0; i < numRowInds; ++i) {
    mMatrix->getLocalRowView(i, localColInds, vals);
    numEntries = localColInds.size();

    for (size_t j = 0; j < numEntries; ++j)
      globColInds[j] = colInds[localColInds[j]];

    this->insertEpetraValues(rowInds[i], numEntries,
        &globColInds[0], &vals[0], res);
  }
  res->FillComplete(*colMap, *rowMap);

  return res;
}
#endif

template class MxCrsMatrix<double>;
template class MxCrsMatrix<MxComplex>;
