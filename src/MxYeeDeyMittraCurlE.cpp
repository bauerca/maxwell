
#include "MxYeeDeyMittraCurlE.h"
#include "MxTypes.h"
#include "MxDimVector.hpp"
#include "MxGridFieldIter.hpp"



template<size_t DIM, typename Scalar>
MxYeeDeyMittraCurlE<DIM, Scalar>::MxYeeDeyMittraCurlE(
RCP<MxEMSim<DIM> > theSim) : 
//MxCrsMatrix<Scalar>(theSim->getField("bfield").getMap(),
//  theSim->getField("efield").getMap()),
MxCrsMatrix<Scalar>(theSim->getField("bfield").getMap()),
efield(&theSim->getField("efield")), 
bfield(&theSim->getField("bfield")),
grid(&theSim->getGrid()) {
  if (DIM == 2)
    setMatrix2d();
  else
    setMatrix3d();
}


template<size_t DIM, typename Scalar>
void MxYeeDeyMittraCurlE<DIM, Scalar>::setMatrix2d() {
  size_t numBComps = bfield->getNumComps();

  MxDimVector<double, DIM> cellSize(grid->getCellSize());
  Scalar idx, idy;
  MxUtil::convertScalar(1 / cellSize[0], idx);
  MxUtil::convertScalar(1 / cellSize[1], idy);

  MxGridFieldIter<DIM> fieldIter(bfield);

  MxIndex row;
  MxIndex cols[4];
  Scalar val, factor;
  //double vals[4] = {1, -1, -1, 1};
  Scalar vals[4] = {idx, -idx, -idy, idy};
  MxComplex factors[4];

  size_t comp0, comp1;
  MxDimVector<int, DIM> cell;

  switch (numBComps) {
  case 1:
    //matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, bfield->getMap(), 4));
    // TE: Bz, Ex, Ey
    for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
      cell = fieldIter.getCell();
      row = bfield->globCompIndx(0, cell);

      cell[0]++;
      cols[0] = efield->globCompIndx(1, cell);
      factors[0] = efield->getCompFactor(1, cell);
      cell[0]--;
      cols[1] = efield->globCompIndx(1, cell);
      factors[1] = efield->getCompFactor(1, cell);
      cell[1]++;
      cols[2] = efield->globCompIndx(0, cell);
      factors[2] = efield->getCompFactor(0, cell);
      cell[1]--;
      cols[3] = efield->globCompIndx(0, cell);
      factors[3] = efield->getCompFactor(0, cell);

      for (size_t i = 0; i < 4; i++) {
        if (factors[i] != Teuchos::ScalarTraits<MxComplex>::zero()) {
          MxUtil::convertScalar(factors[i], factor);
          val = factor * vals[i];
          MxCrsMatrix<Scalar>::insertRowValues(row, 1, &cols[i], &val);
        }
      }
    }
    break;
  case 2:
    //matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, bfield->getMap(), 2));
    // TM: Ez, Bx, By
    for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
      cell = fieldIter.getCell();

      comp0 = fieldIter.getComp(); // Bx : 0; By : 1
      comp1 = (comp0 + 1) % 2;

      if (bfield->getCompFactor(comp0, cell) == Teuchos::ScalarTraits<MxComplex>::zero()) {
        continue;
      }

      row = bfield->globCompIndx(comp0, cell);

      cell[comp1]++;
      cols[0] = efield->globCompIndx(0, cell);
      factors[0] = efield->getCompFactor(0, cell);
      cell[comp1]--;
      cols[1] = efield->globCompIndx(0, cell);
      factors[1] = efield->getCompFactor(0, cell);

      for (size_t i = 0; i < 2; i++) {
        if (factors[i] != Teuchos::ScalarTraits<MxComplex>::zero()) {
          MxUtil::convertScalar(factors[i], factor);
          val = factor * vals[i + 2 * comp1];
          MxCrsMatrix<Scalar>::insertRowValues(row, 1, &cols[i], &val);
        }
      }
    }
    MxCrsMatrix<Scalar>::scale(-Teuchos::ScalarTraits<Scalar>::one());
    break;
  }
  MxCrsMatrix<Scalar>::fillComplete(efield->getMap(), bfield->getMap());
}


// (curl E)_x = dEz/dy - dEy/dz
// (curl E)_y = dEx/dz - dEz/dx
// (curl E)_z = dEy/dx - dEx/dy = (Ey^+ - Ey^-)/dx + (
template<size_t DIM, typename Scalar>
void MxYeeDeyMittraCurlE<DIM, Scalar>::setMatrix3d() {
  //matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, bfield->getMap(), 4));

  MxDimVector<double, DIM> cellSize(grid->getCellSize());
  Scalar idx, idy, idz;
  MxUtil::convertScalar(1 / cellSize[0], idx);
  MxUtil::convertScalar(1 / cellSize[1], idy);
  MxUtil::convertScalar(1 / cellSize[2], idz);

  MxIndex row;
  MxIndex cols[4];
  Scalar val, factor;
  //double vals[4] = {1, -1, -1, 1};
  Scalar vals[3][4] = {{idy, -idy, -idz, idz}, 
                       {idz, -idz, -idx, idx}, 
                       {idx, -idx, -idy, idy}};
  MxComplex factors[4];

  size_t comp0, comp1, comp2;
  MxDimVector<int, DIM> cell;

  MxGridFieldIter<DIM> fieldIter(bfield);

  for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
    cell = fieldIter.getCell();
    //std::cout << cell;

    comp0 = fieldIter.getComp();
    comp1 = (comp0 + 1) % 3;
    comp2 = (comp0 + 2) % 3;

    if (bfield->getCompFactor(comp0, cell) == Teuchos::ScalarTraits<MxComplex>::zero()) {
      continue;
    }

    row = bfield->globCompIndx(comp0, cell);

    cell[comp1]++;
    cols[0] = efield->globCompIndx(comp2, cell);
    factors[0] = efield->getCompFactor(comp2, cell);
    cell[comp1]--;
    cols[1] = efield->globCompIndx(comp2, cell);
    factors[1] = efield->getCompFactor(comp2, cell);
    cell[comp2]++;
    cols[2] = efield->globCompIndx(comp1, cell);
    factors[2] = efield->getCompFactor(comp1, cell);
    cell[comp2]--;
    cols[3] = efield->globCompIndx(comp1, cell);
    factors[3] = efield->getCompFactor(comp1, cell);

    for (size_t i = 0; i < 4; i++) {
      if (factors[i] != Teuchos::ScalarTraits<MxComplex>::zero()) {
        MxUtil::convertScalar(factors[i], factor);
        val = factor * vals[comp0][i];
        MxCrsMatrix<Scalar>::insertRowValues(row, 1, &cols[i], &val);
      }
    }
    //std::cout << "\n";

  }
  //this->FillComplete(efield->getMap(), bfield->getMap());
  MxCrsMatrix<Scalar>::fillComplete(efield->getMap(), bfield->getMap());
}

template class MxYeeDeyMittraCurlE<1, double>;
template class MxYeeDeyMittraCurlE<2, double>;
template class MxYeeDeyMittraCurlE<3, double>;
template class MxYeeDeyMittraCurlE<1, MxComplex>;
template class MxYeeDeyMittraCurlE<2, MxComplex>;
template class MxYeeDeyMittraCurlE<3, MxComplex>;
