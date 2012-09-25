
#include "MxYeeDeyMittraFracs.h"
#include "MxDimVector.hpp"
#include "MxGridFieldIter.hpp"

#include <cstdlib>



template<size_t DIM, typename Scalar>
MxYeeDeyMittraFracs<DIM, Scalar>::MxYeeDeyMittraFracs(
const MxGridField<DIM> * aGridField, 
const MxGrid<DIM> * aGrid, bool getInverse, double aMinFrac, bool randomize,
double randomScale) : 
MxCrsMatrix<Scalar>(aGridField->getMap()),
mRandomize(randomize),
mRandomScale(randomScale),
minFrac(aMinFrac),
inverse(getInverse),
field(aGridField),
grid(aGrid) {
  if (inverse)
    setMatrixInverse();
  else
    setMatrix();
}


template<size_t DIM, typename Scalar>
void MxYeeDeyMittraFracs<DIM, Scalar>::setMatrix() {
  //matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, bfield->getMap(), 1, true));
  //

  MxGridFieldIter<DIM> fieldIter(field);

  double minNonzeroFrac = 1.0;
  int numZeroFracs = 0;

  MxIndex ind;
  size_t comp;
  double val;
  Scalar mval;
  MxDimVector<int, DIM> cell;
  for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
    ind = fieldIter.getGlobCompIndx();
    comp = fieldIter.getComp();
    cell = fieldIter.getCell();
    val = field->getCompFrac(comp, cell, "pec");
    //std::cout << val << "\n";
    if (val != 0) {
      if (val < minFrac) val = minFrac;
      if (val < minNonzeroFrac) minNonzeroFrac = val;
    }
    else
      numZeroFracs++;

    MxUtil::convertScalar(val, mval);
    MxCrsMatrix<Scalar>::insertRowValues(ind, 1, &ind, &mval);
  }

  std::cout << "  min nonzero frac is: " << minNonzeroFrac << "\n";
  std::cout << "  number of zeros is: " << numZeroFracs << "\n";

  //matrix->FillComplete(bfield->getMap(), bfield->getMap());
  MxCrsMatrix<Scalar>::fillComplete(field->getMap(), field->getMap());
}

template<size_t DIM, typename Scalar>
void MxYeeDeyMittraFracs<DIM, Scalar>::setMatrixInverse() {
  //matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, bfield->getMap(), 1, true));
  //

  MxGridFieldIter<DIM> fieldIter(field);

  double minNonzeroFrac = 1.0;
  int numZeroFracs = 0;

  MxIndex ind;
  size_t comp;
  double val;
  Scalar mval;
  MxDimVector<int, DIM> cell;

  srand(time(NULL) + this->getRangeMap()->getComm()->myPID());
  double rnd, df, nval;

  for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
    ind = fieldIter.getGlobCompIndx();
    comp = fieldIter.getComp();
    cell = fieldIter.getCell();
    val = field->getCompFrac(comp, cell, "pec");

    // pertub testing
    if (mRandomize) {
      if (val != 0 && val != 1) {
        rnd = (2 * double(rand() - RAND_MAX / 2)) / double(RAND_MAX);
        df = mRandomScale * rnd;
        //nval = val + df;
        nval = val*pow(10.0, df);
        if (nval < 0) val = 0;
        else if (nval > 1) val = 1;
        else {
          //std::cout << "original=" << val << " new=" << nval << "\n";
          val = nval;
        }
      }
    }

    if (val == 0) 
      numZeroFracs++;
    else if (val < minFrac) {
      val = 1.0 / minFrac;
      if (val < minNonzeroFrac) minNonzeroFrac = minFrac;
    }
    else {
      if (val < minNonzeroFrac) minNonzeroFrac = val;
      val = 1.0 / val;
    }

    MxUtil::convertScalar(val, mval);
    MxCrsMatrix<Scalar>::insertRowValues(ind, 1, &ind, &mval);
  }

  std::cout << "  min nonzero frac is: " << minNonzeroFrac << "\n";
  std::cout << "  number of zeros is: " << numZeroFracs << "\n";

  //matrix->FillComplete(bfield->getMap(), bfield->getMap());
  MxCrsMatrix<Scalar>::fillComplete(field->getMap(), field->getMap());
}

template class MxYeeDeyMittraFracs<1, double>;
template class MxYeeDeyMittraFracs<2, double>;
template class MxYeeDeyMittraFracs<3, double>;
template class MxYeeDeyMittraFracs<1, MxComplex>;
template class MxYeeDeyMittraFracs<2, MxComplex>;
template class MxYeeDeyMittraFracs<3, MxComplex>;
