#ifndef MX_MULTI_VECTOR
#define MX_MULTI_VECTOR

#include "MxTypes.h"

#include "Epetra_MultiVector.h"
#include "Tpetra_MultiVector.hpp"

class MxMap;
template<typename Scalar> class MxCrsMatrix;
template<typename Scalar> class MxVector;

template<typename Scalar>
class MxMultiVector {
  public:
    MxMultiVector(RCP<MxMap> map, size_t numVecs);

    // ctor
    MxMultiVector(MxMultiVector<Scalar> const & mv);

    MxMultiVector(MxMultiVector<Scalar> const & mv,
        std::vector<size_t> const & vecInds, bool deepcopy);

    size_t getNumVecs() const {return mNumVecs;}

    /**
     *  Indices are always global here 
     */
    void replaceGlobalValue(MxIndex row, size_t vec, Scalar value);

    void replaceLocalValue(MxIndex row, size_t vec, Scalar value);

    //void scaleLocalValue(MxIndex row, size_t vec, Scalar value);

    void conj();

    void scale(Scalar val);

    void scale(std::vector<Scalar> const & vals);

    virtual void random() {mVec->Random();}

    void norm2(std::vector<double> & norms) const;

    void normalize();

    /**
    Sets all values in multivector to val
    */
    void set(Scalar val);

    /**
    Computes mv^{\dagger} \cdot *this
    */
    void dot(MxMultiVector<Scalar> const & mv, std::vector<Scalar> & res) const;

    void update(Scalar scalarA, MxMultiVector<Scalar> const & mvA, Scalar scalarThis);

    void setVector(size_t vecIndex, RCP<MxVector<Scalar> const> vec, bool copy);

    RCP<MxVector<Scalar> const> getVector(size_t vecIndex, bool copy) const;

    RCP<MxVector<Scalar> > getVectorNonConst(size_t vecIndex, bool copy);

    Scalar operator()(MxIndex localIndex, size_t vec) const;

    MxMultiVector<Scalar> & operator=(MxMultiVector<Scalar> const & mv);


#if 0
    RCP<Tpetra::MultiVector<Scalar, MxIndex> > getTpetraMultiVector() const {
      return mVec;
    }

    RCP<Epetra_MultiVector> getEpetraMultiVector();
#endif

    RCP<MxMap> getMap() const {return mMap;}

    MxIndex getLocalLength() const;

// advanced usage!
#if LINALG_BASE == TPETRA
#elif LINALG_BASE == EPETRA
    RCP<Epetra_MultiVector> getRawMV() {return mVec;}

    RCP<const Epetra_MultiVector> getRawMV() const {return mVec;}

    void setRawMV(RCP<Epetra_MultiVector> const & mv) {mVec = mv;}
#endif

    friend class MxCrsMatrix<Scalar>;

  private:
    size_t mNumVecs;

    RCP<MxMap> mMap;

#if LINALG_BASE == TPETRA
    RCP<Tpetra::MultiVector<Scalar, MxIndex> > mVec;
#elif LINALG_BASE == EPETRA
    RCP<Epetra_MultiVector> mVec;

    RCP<MxCrsMatrix<Scalar> > mEye;
#endif

    //RCP<Epetra_MultiVector> mEpetraVec;

};


#if LINALG_BASE == TPETRA
template<typename Scalar>
inline
void MxMultiVector<Scalar>::replaceLocalValue(MxIndex row, size_t vec, Scalar value) {
  mVec->replaceLocalValue(row, vec, value);
}
#elif LINALG_BASE == EPETRA
template<>
inline
void MxMultiVector<double>::replaceLocalValue(MxIndex row, size_t vec, double value) {
  mVec->ReplaceMyValue(int(row), int(vec), value);
}

template<>
inline
void MxMultiVector<MxComplex>::replaceLocalValue(MxIndex row, size_t vec, MxComplex value) {
  mVec->ReplaceMyValue(int(2*row), int(vec), value.real());
  mVec->ReplaceMyValue(int(2*row+1), int(vec), value.imag());
}
#endif



#if LINALG_BASE == TPETRA
template<typename Scalar>
inline
void MxMultiVector<Scalar>::replaceGlobalValue(MxIndex row, size_t vec, Scalar value) {
  mVec->replaceGlobalValue(row, vec, value);
}
#elif LINALG_BASE == EPETRA
template<>
inline
void MxMultiVector<double>::replaceGlobalValue(MxIndex row, size_t vec, double value) {
  mVec->ReplaceGlobalValue(int(row), int(vec), value);
}

template<>
inline
void MxMultiVector<MxComplex>::replaceGlobalValue(MxIndex row, size_t vec, MxComplex value) {
  mVec->ReplaceGlobalValue(int(2*row), int(vec), value.real());
  mVec->ReplaceGlobalValue(int(2*row+1), int(vec), value.imag());
}
#endif





#if LINALG_BASE == TPETRA
template<typename Scalar>
inline
Scalar MxMultiVector<Scalar>::operator()(MxIndex localIndex, size_t vec) const {
  return mVec->getData(vec)[localIndex];
}
#elif LINALG_BASE == EPETRA
template<>
inline
double MxMultiVector<double>::operator()(MxIndex row, size_t vec) const {
  return mVec->operator[](vec)[row];
}

template<>
inline
MxComplex MxMultiVector<MxComplex>::operator()(MxIndex row, size_t vec) const {
  double re, im;
  //std::cout << (2*row) << ", my lid? " << mVec->Map().MyLID(2*row) << "\n";
  //std::cout << (2*row+1) << ", my lid? " << mVec->Map().MyLID(2*row+1) << "\n";
  re = mVec->operator[](vec)[2*row];
  im = mVec->operator[](vec)[2*row+1];
  return MxComplex(re, im);
}
#endif


#if LINALG_BASE == TPETRA
#elif LINALG_BASE == EPETRA
template<>
inline
MxIndex MxMultiVector<double>::getLocalLength() const {
  return mVec->MyLength();
}

template<>
inline
MxIndex MxMultiVector<MxComplex>::getLocalLength() const {
  return mVec->MyLength() / 2;
}
#endif


#endif
