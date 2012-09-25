#ifndef MX_CRS_MATRIX
#define MX_CRS_MATRIX

#include "MxTypes.h"
#include "MxVector.hpp"

#include "Tpetra_CrsMatrix.hpp"

class Epetra_CrsMatrix;

template<typename Scalar>
class MxCrsMatrix {
  public:
    MxCrsMatrix(RCP<MxMap> rowMap);
    
    MxCrsMatrix(RCP<MxMap> rowMap, RCP<MxMap> colMap);

    /**
    Create a diagonal matrix from a scalar
    */
    MxCrsMatrix(RCP<MxMap> rowMap, Scalar diag);

    /**
    Create a diagonal matrix from a vector
    */
    MxCrsMatrix(RCP<MxVector<Scalar> > diag);

    ~MxCrsMatrix() {};

    /**
     *  Indices are always global here 
     */
    void insertRowValues(MxIndex row, std::vector<MxIndex> cols, std::vector<Scalar> vals);

    void insertRowValues(MxIndex row, size_t numEntries, MxIndex const * cols, Scalar const * vals);

    void scale(Scalar val);

    void leftScale(RCP<MxVector<Scalar> > vec);

    void getLocalRow(MxIndex localRowIndex,
      std::vector<MxIndex> & localColIndices, std::vector<Scalar> & vals) const;

    RCP<MxVector<Scalar> > getRowSums(bool inverse) const;

    void apply(MxMultiVector<Scalar> const & x, MxMultiVector<Scalar> & y) const;

#if LINALG_BASE == TPETRA
#elif LINALG_BASE == EPETRA
    RCP<Epetra_CrsMatrix> getRawMatrix() {return mMatrix;}

    RCP<Epetra_CrsMatrix const> getRawMatrix() const {return mMatrix;}
#endif

    RCP<Tpetra::CrsMatrix<Scalar, MxIndex> > getTpetraCrsMatrix() const;

    void save(std::string name) const;

    void purgeZeros();

    /**
     *  FillComplete should have been called first
     */
    RCP<Epetra_CrsMatrix> getEpetraCrsMatrix() const;

    void fillComplete(RCP<MxMap> domainMap, RCP<MxMap> rangeMap);

    bool isFilled() const {return mFilled;}

    RCP<MxMap> getDomainMap() const {return mColMap;}

    RCP<MxMap> getRangeMap() const {return mRowMap;}

    static void matrixMatrixMultiply(
      MxCrsMatrix<Scalar> const & a, bool transA,
      MxCrsMatrix<Scalar> const & b, bool transB,
      MxCrsMatrix<Scalar> & result, bool fc);

    static void matrixMatrixAdd(
      MxCrsMatrix<Scalar> const & a, bool transA, Scalar scalarA,
      MxCrsMatrix<Scalar> const & b, bool transB, Scalar scalarB,
      MxCrsMatrix<Scalar> & result, bool fc);

  private:
    bool mFilled;

    RCP<MxMap> mRowMap, mColMap;

    RCP<Epetra_Map> mEpetraRowMap, mEpetraColMap;

#if LINALG_BASE == TPETRA
    RCP<Tpetra::CrsMatrix<Scalar, MxIndex> > mMatrix;
#elif LINALG_BASE == EPETRA
    RCP<Epetra_CrsMatrix> mMatrix;
#endif

    void initBaseMatrix();

    void insertEpetraValues(MxIndex row, size_t numEntries,
        MxIndex const * inds, Scalar const * vals,
        RCP<Epetra_CrsMatrix> matrix) const;

};

#endif
