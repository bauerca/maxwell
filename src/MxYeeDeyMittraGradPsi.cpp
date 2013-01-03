
#include "MxYeeDeyMittraGradPsi.h"
#include "MxTypes.h"
#include "MxDimVector.hpp"
#include "MxGridFieldIter.hpp"



template<size_t DIM, typename Scalar>
MxYeeDeyMittraGradPsi<DIM, Scalar>::MxYeeDeyMittraGradPsi(
RCP<MxEMSim<DIM> > theSim) :
//MxCrsMatrix<Scalar>(theSim->getField("bfield").getMap(),
//  theSim->getField("psifield").getMap()),
MxCrsMatrix<Scalar>(theSim->getField("bfield").getMap()),
psifield(&theSim->getField("psifield")),
bfield(&theSim->getField("bfield")),
grid(&theSim->getGrid()) {
  setMatrix();
}


template<size_t DIM, typename Scalar>
void MxYeeDeyMittraGradPsi<DIM, Scalar>::setMatrix() {

  size_t numBComps = bfield->getNumComps();

  if (numBComps == 1) {
    std::cout << "MxYeeDeyMittraGradPsi::setMatrix(): 2D TE simulation, gradient of Psi "
              << "is always zero. This matrix is zero.\n";
    //matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, psifield->getMap(), 0));
  }
  else {
    //matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, psifield->getMap(), 4));

    MxDimVector<double, DIM> invCellSize(1.0 / grid->getCellSize());

    MxDimVector<MxDimVector<double, 2>, DIM> vals;
    for (size_t i = 0; i < DIM; ++i) {
      MxUtil::convertScalar(invCellSize[i], vals[i][0]);
      MxUtil::convertScalar(-invCellSize[i], vals[i][1]);
      //vals[i][0] = invCellSize[i];
      //vals[i][1] = -invCellSize[i];
    }

    MxGridFieldIter<DIM> fieldIter(bfield);

    MxIndex row;
    size_t comp;
    Scalar val, factor;
    MxIndex cols[2];
    MxComplex factors[2];

    MxDimVector<int, DIM> cell;

    for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
      cell = fieldIter.getCell();
      comp = fieldIter.getComp();

      if (bfield->getCompFactor(comp, cell) == ScalarTraits<MxComplex>::zero()) {
        continue;
      }

      row = bfield->globCompIndx(comp, cell);

      cols[0] = psifield->globCompIndx(0, cell);
      factors[0] = psifield->getCompFactor(0, cell);
      cell[comp]--;
      cols[1] = psifield->globCompIndx(0, cell);
      factors[1] = psifield->getCompFactor(0, cell);

      for (size_t i = 0; i < 2; i++) {
        //std::cout << factors[i] << ", " << cols[i] << "\n";
        if (factors[i] != ScalarTraits<MxComplex>::zero()) {
          MxUtil::convertScalar(factors[i], factor);
          val = factor * vals[comp][i];
          //std::cout << "val = " << val << "\n";
          MxCrsMatrix<Scalar>::insertRowValues(row, 1, &cols[i], &val);
        }
      }
    }
    MxCrsMatrix<Scalar>::fillComplete(psifield->getMap(), bfield->getMap());

  }
}


template class MxYeeDeyMittraGradPsi<1, double>;
template class MxYeeDeyMittraGradPsi<2, double>;
template class MxYeeDeyMittraGradPsi<3, double>;
template class MxYeeDeyMittraGradPsi<1, MxComplex>;
template class MxYeeDeyMittraGradPsi<2, MxComplex>;
template class MxYeeDeyMittraGradPsi<3, MxComplex>;
