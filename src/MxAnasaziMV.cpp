#include "MxAnasaziMV.hpp"

#include "MxComm.hpp"

#include "Epetra_LocalMap.h"

template<>
void MxAnasaziMV<double>::MvTimesMatAddMv(double alpha,
    const Anasazi::MultiVec<double> & A,
    const Teuchos::SerialDenseMatrix<int,double> & B, double beta) {
  //std::cout << "MvTimesMatAddMv called\n";
  Epetra_LocalMap LocalMap(
      int(B.numRows()), 0, *this->getMap()->getComm()->getEpetraComm());
  Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());

  MxAnasaziMV<double> *A_vec =
      dynamic_cast<MxAnasaziMV<double> *>(
          &const_cast<Anasazi::MultiVec<double> &>(A));

  if (A_vec == NULL) {
    std::cout << "MxAnasaziMV<double>::MvTimesMatAddMv cast failed\n";
    exit(EXIT_FAILURE);
  }

  int res = this->getRawMV()->Multiply(
      'N', 'N', alpha, *A_vec->getRawMV(), B_Pvec, beta);

  if (res != 0) {
    std::cout << "MxAnasaziMV<double>::MvTimesMatAddMv call to \
      Epetra_MultiVector::Multiply failed.\n";
    exit(EXIT_FAILURE);
  }
}


/*
 * Performs   this = alpha*A*B + beta*this
 */
template<>
void MxAnasaziMV<MxComplex>::MvTimesMatAddMv(MxComplex alpha,
    const Anasazi::MultiVec<MxComplex> & A,
    const Teuchos::SerialDenseMatrix<int, MxComplex> & B,
    MxComplex beta) {

  //std::cout << "MvTimesMatAddMv called\n";
  size_t numRows = B.numRows();
  size_t numCols = B.numCols();
  Teuchos::SerialDenseMatrix<int, double> reB(numRows, numCols);
  Teuchos::SerialDenseMatrix<int, double> imB(numRows, numCols);
  for (size_t i = 0; i < numRows; ++i) {
    for (size_t j = 0; j < numCols; ++j) {
      reB(i, j) = B(i, j).real();
      imB(i, j) = B(i, j).imag();
    }
  }
  Epetra_LocalMap lMap(int(numRows), 0, *getMap()->getComm()->getEpetraComm());
  Epetra_MultiVector reBView(View, lMap, reB.values(), reB.stride(), numCols);
  Epetra_MultiVector imBView(View, lMap, imB.values(), imB.stride(), numCols);

  //std::cout << alpha << ", " << beta << "\n";
  //std::cout << B;
  //std::cout << reBView << imBView;

  MxAnasaziMV<MxComplex> * castedA =
      dynamic_cast<MxAnasaziMV<MxComplex> *>(
          &const_cast<Anasazi::MultiVec<MxComplex> &>(A));

  MxMultiVector<MxComplex> AxBRe(*this), AxBIm(*this);
  AxBRe.getRawMV()->Multiply(
      'N', 'N', 1.0, *castedA->getRawMV(), reBView, 0.0);
  AxBIm.getRawMV()->Multiply(
      'N', 'N', 1.0, *castedA->getRawMV(), imBView, 0.0);
  AxBIm.scale(MxUtil::i<MxComplex>());
  //std::cout << "A*Bre\n";
  //std::cout << *AxBRe.getRawMV();
  //std::cout << "A*Bim\n";
  //std::cout << *AxBIm.getRawMV();
  
  MxMultiVector<MxComplex> & AxB = AxBRe;
  AxB.update(alpha, AxBIm, alpha);
  
  //std::cout << "alpha*A*B\n";
  //std::cout << *AxB.getRawMV();

  this->update(MxComplex(1.0, 0.0), AxB, beta);
}

template<typename Scalar>
void MxAnasaziMV<Scalar>::MvAddMv(Scalar alpha,
    const Anasazi::MultiVec<Scalar> & A, Scalar beta,
    const Anasazi::MultiVec<Scalar> & B)
{
  //std::cout << "MvAddMv called\n";
  MxAnasaziMV<Scalar> * ptrA =        
    dynamic_cast<MxAnasaziMV<Scalar> *>(
        &const_cast<Anasazi::MultiVec<Scalar> &>(A));
  MxAnasaziMV<Scalar> * ptrB =        
    dynamic_cast<MxAnasaziMV<Scalar> *>(
        &const_cast<Anasazi::MultiVec<Scalar> &>(B));

  //std::cout << "this\n" << *this->getRawMV();
  //std::cout << "A\n" << *ptrA->getRawMV();
  //std::cout << "B\n" << *ptrB->getRawMV();

  //std::cout << "before\n" << *this->getRawMV();
  // Assumption that self-assignment of underlying multivector is okay. This
  // is up to the Trilinos people. We trust them!
  *this->getRawMV() = *ptrA->getRawMV();
  //std::cout << "after\n" << *this->getRawMV();
  this->update(beta, *ptrB, alpha);
}

template<>
void MxAnasaziMV<double>::MvTransMv(double alpha,
    const Anasazi::MultiVec<double> & A,
    Teuchos::SerialDenseMatrix<int, double> & B) const {
  //std::cout << "MvTransMv called\n";
  //std::cout << B.numRows() << "x" << B.numCols() << "\n";
  Epetra_LocalMap LocalMap(
      int(B.numRows()), 0, *this->getMap()->getComm()->getEpetraComm());
  Epetra_MultiVector mvB(View, LocalMap, B.values(), B.stride(), B.numCols());

  MxAnasaziMV<double> *A_vec =
      dynamic_cast<MxAnasaziMV<double> *>(
          &const_cast<Anasazi::MultiVec<double> &>(A));

  if (A_vec == NULL) {
    std::cout << "MxAnasaziMV<double>::MvTransMv cast failed\n";
    exit(EXIT_FAILURE);
  }

  //std::cout << "B: " << mvB.GlobalLength() << "x" << mvB.NumVectors() << "\n";
  //std::cout << "A: " <<
  //    A_vec->getRawMV()->GlobalLength() << "x" <<
  //    A_vec->getRawMV()->NumVectors() << "\n";
  //std::cout << "this: " <<
  //    this->getRawMV()->GlobalLength() << "x" <<
  //    this->getRawMV()->NumVectors() << "\n";
  int res = mvB.Multiply(
      'T', 'N', alpha, *A_vec->getRawMV(), *this->getRawMV(), 0.0);

  if (res != 0) {
    std::cout << "MxAnasaziMV<double>::MvTransMv call to \
      Epetra_MultiVector::Multiply failed.\n";
    std::cout << "Error code: " << res << "\n";
    exit(EXIT_FAILURE);
  }
}


template<>
void MxAnasaziMV<MxComplex>::MvTransMv(MxComplex alpha,
    const Anasazi::MultiVec<MxComplex> & A,
    Teuchos::SerialDenseMatrix<int, MxComplex> & B) const {

  //std::cout << "alpha: " << alpha << "\n";
  //std::cout << "this: " << this->GetVecLength() << "x" << this->GetNumberVecs() << "\n";
  //std::cout << "A: " << A.GetVecLength() << "x" << A.GetNumberVecs() << "\n";
  //std::cout << "B: " << B.numRows() << "x" << B.numCols() << "\n";

  size_t numRows = B.numRows();
  size_t numCols = B.numCols();
  Teuchos::SerialDenseMatrix<int, double> reB(numRows, numCols);
  Teuchos::SerialDenseMatrix<int, double> imB(numRows, numCols);
  for (size_t i = 0; i < numRows; ++i) {
    for (size_t j = 0; j < numCols; ++j) {
      reB(i, j) = B(i, j).real();
      imB(i, j) = B(i, j).imag();
    }
  }
  Epetra_LocalMap lMap(int(numRows), 0, *getMap()->getComm()->getEpetraComm());
  Epetra_MultiVector reBView(View, lMap, reB.values(), reB.stride(), numCols);
  Epetra_MultiVector imBView(View, lMap, imB.values(), imB.stride(), numCols);

  MxAnasaziMV<MxComplex> * castedA =
      dynamic_cast<MxAnasaziMV<MxComplex> *>(
          &const_cast<Anasazi::MultiVec<MxComplex> &>(A));
  // Do this for adjoint
  MxMultiVector<MxComplex> iA(*castedA);
  iA.scale(MxUtil::i<MxComplex>());

  reBView.Multiply('T', 'N', 1.0, *castedA->getRawMV(), *this->getRawMV(), 0.0);
  imBView.Multiply('T', 'N', 1.0, *iA.getRawMV(), *this->getRawMV(), 0.0);

  // Do this for transpose
  //MxMultiVector<MxComplex> imA(*castedA), reA(*castedA);
  //reA.conj();
  //imA.conj();
  //imA.scale(MxUtil::i<MxComplex>());

  //reBView.Multiply('T', 'N', 1.0, *reA.getRawMV(), *this->getRawMV(), 0.0);
  //imBView.Multiply('T', 'N', 1.0, *imA.getRawMV(), *this->getRawMV(), 0.0);

  for (size_t i = 0; i < numRows; ++i)
    for (size_t j = 0; j < numCols; ++j)
      B(i, j) = alpha * MxComplex(reB(i, j), imB(i, j));
}


template<typename Scalar>
void MxAnasaziMV<Scalar>::SetBlock(const Anasazi::MultiVec<Scalar> & A,
    const std::vector<int> & index) {
  //std::cout << "SetBlock called\n";
  int n = this->getRawMV()->Map().NumMyElements();
  int nvecs = A.GetNumberVecs();
  Epetra_MultiVector * ptrThis = this->getRawMV().get();
  Epetra_MultiVector * ptrA = 
      dynamic_cast<MxAnasaziMV<Scalar> &>(
          const_cast<Anasazi::MultiVec<Scalar> &>(A)).getRawMV().get();
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < nvecs; ++j)
      ptrThis->ReplaceMyValue(i, index[j], (*ptrA)[j][i]);
}




template class MxAnasaziMV<double>;
template class MxAnasaziMV<MxComplex>;
