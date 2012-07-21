
#include "MxYeeDeyMittraCurlB.h"
#include "MxTypes.h"
#include "MxDimVector.hpp"
#include "MxGridFieldIter.hpp"



template<size_t DIM, typename Scalar>
MxYeeDeyMittraCurlB<DIM, Scalar>::MxYeeDeyMittraCurlB(
RCP<MxEMSim<DIM> > theSim) :
//MxCrsMatrix<Scalar>(theSim->getField("efield").getMap(),
//  theSim->getField("bfield").getMap()),
MxCrsMatrix<Scalar>(theSim->getField("efield").getMap()),
efield(&theSim->getField("efield")), 
bfield(&theSim->getField("bfield")), 
grid(&theSim->getGrid()) {
  if (DIM == 2)
    setMatrix2d();
  else
    setMatrix3d();
}


template<size_t DIM, typename Scalar>
void MxYeeDeyMittraCurlB<DIM, Scalar>::setMatrix2d() {
  size_t numEComps = efield->getNumComps();

  MxDimVector<double, DIM> cellSize(grid->getCellSize());
  Scalar idx, idy;
  MxUtil::convertScalar(1 / cellSize[0], idx);
  MxUtil::convertScalar(1 / cellSize[1], idy);

  MxGridFieldIter<DIM> fieldIter(efield);

  MxIndex row;
  Scalar val, factor;

  MxIndex cols[4];
  //double vals[4] = {1, -1, -1, 1};
  Scalar vals[4] = {idx, -idx, -idy, idy}; //implicit cast to complex
  MxComplex factors[4];

  size_t comp0, comp1;
  MxDimVector<int, DIM> cell;

  switch (numEComps) {
  case 1:
    //matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, bfield->getMap(), 4));
    // TM: Ez, Bx, By
    for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
      cell = fieldIter.getCell();
      row = efield->globCompIndx(0, cell);

      cols[0] = bfield->globCompIndx(1, cell);
      factors[0] = bfield->getCompFactor(1, cell);
      cell[0]--;
      cols[1] = bfield->globCompIndx(1, cell);
      factors[1] = bfield->getCompFactor(1, cell);
      cell[0]++;
      cols[2] = bfield->globCompIndx(0, cell);
      factors[2] = bfield->getCompFactor(0, cell);
      cell[1]--;
      cols[3] = bfield->globCompIndx(0, cell);
      factors[3] = bfield->getCompFactor(0, cell);

      for (size_t i = 0; i < 4; i++) {
        if (factors[i] != ScalarTraits<MxComplex>::zero()) {
          MxUtil::convertScalar(factors[i], factor);
          val = factor * vals[i];
          MxCrsMatrix<Scalar>::insertRowValues(row, 1, &cols[i], &val);
        }
      }
    }
    break;
  case 2:
    //matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, bfield->getMap(), 2));
    // TE: Bz, Ex, Ey
    for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
      cell = fieldIter.getCell();

      comp0 = fieldIter.getComp(); // Bx : 0; By : 1
      comp1 = (comp0 + 1) % 2;

      row = efield->globCompIndx(comp0, cell);

      cols[0] = bfield->globCompIndx(0, cell);
      factors[0] = bfield->getCompFactor(0, cell);
      cell[comp1]--;
      cols[1] = bfield->globCompIndx(0, cell);
      factors[1] = bfield->getCompFactor(0, cell);

      for (size_t i = 0; i < 2; i++) {
        if (factors[i] != ScalarTraits<MxComplex>::zero()) {
          MxUtil::convertScalar(factors[i], factor);
          val = factor * vals[i + 2 * comp1];
          MxCrsMatrix<Scalar>::insertRowValues(row, 1, &cols[i], &val);
        }
      }
    }
    MxCrsMatrix<Scalar>::scale(-ScalarTraits<Scalar>::one());
    break;
  }

  MxCrsMatrix<Scalar>::fillComplete(bfield->getMap(), efield->getMap());

}

// (curl E)_x = dEz/dy - dEy/dz
// (curl E)_y = dEx/dz - dEz/dx
// (curl E)_z = dEy/dx - dEx/dy = (Ey^+ - Ey^-)/dx + (
template<size_t DIM, typename Scalar>
void MxYeeDeyMittraCurlB<DIM, Scalar>::setMatrix3d() {
  //matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, bfield->getMap(), 4));

  MxDimVector<double, DIM> cellSize(grid->getCellSize());
  Scalar idx, idy, idz;
  MxUtil::convertScalar(1.0 / cellSize[0], idx);
  MxUtil::convertScalar(1.0 / cellSize[1], idy);
  MxUtil::convertScalar(1.0 / cellSize[2], idz);

  MxIndex row;
  Scalar val, factor;
  MxIndex cols[4];
  //double vals[4] = {1, -1, -1, 1};
  Scalar vals[3][4] = {{idy, -idy, -idz, idz}, 
                       {idz, -idz, -idx, idx}, 
                       {idx, -idx, -idy, idy}};
  MxComplex factors[4];

  size_t comp0, comp1, comp2;
  MxDimVector<int, DIM> cell;

  MxGridFieldIter<DIM> fieldIter(efield);

  for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
    cell = fieldIter.getCell();

    comp0 = fieldIter.getComp();
    comp1 = (comp0 + 1) % 3;
    comp2 = (comp0 + 2) % 3;

    row = efield->globCompIndx(comp0, cell);

    cols[0] = bfield->globCompIndx(comp2, cell);
    factors[0] = bfield->getCompFactor(comp2, cell);
    cell[comp1]--;
    cols[1] = bfield->globCompIndx(comp2, cell);
    factors[1] = bfield->getCompFactor(comp2, cell);
    cell[comp1]++;
    cols[2] = bfield->globCompIndx(comp1, cell);
    factors[2] = bfield->getCompFactor(comp1, cell);
    cell[comp2]--;
    cols[3] = bfield->globCompIndx(comp1, cell);
    factors[3] = bfield->getCompFactor(comp1, cell);

    for (size_t i = 0; i < 4; i++) {
      //std::cout << factors[i] << ", ";
      if (factors[i] != ScalarTraits<MxComplex>::zero()) {
        MxUtil::convertScalar(factors[i], factor);
        val = factor * vals[comp0][i];
        MxCrsMatrix<Scalar>::insertRowValues(row, 1, &cols[i], &val);
      }
    }
    //std::cout << "\n";
  }

  MxCrsMatrix<Scalar>::fillComplete(bfield->getMap(), efield->getMap());
}

template class MxYeeDeyMittraCurlB<1, double>;
template class MxYeeDeyMittraCurlB<2, double>;
template class MxYeeDeyMittraCurlB<3, double>;
template class MxYeeDeyMittraCurlB<1, MxComplex>;
template class MxYeeDeyMittraCurlB<2, MxComplex>;
template class MxYeeDeyMittraCurlB<3, MxComplex>;
