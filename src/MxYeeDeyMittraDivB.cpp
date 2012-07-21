
#include "MxYeeDeyMittraDivB.h"
#include "MxTypes.h"
#include "MxDimVector.hpp"
#include "MxGridFieldIter.hpp"



template<size_t DIM, typename Scalar>
MxYeeDeyMittraDivB<DIM, Scalar>::MxYeeDeyMittraDivB(
RCP<MxEMSim<DIM> > theSim) : 
//MxCrsMatrix<Scalar>(theSim->getField("psifield").getMap(),
//  theSim->getField("bfield").getMap()),
MxCrsMatrix<Scalar>(theSim->getField("psifield").getMap()),
bfield(&theSim->getField("bfield")),
psifield(&theSim->getField("psifield")),
grid(&theSim->getGrid()) {
  if (bfield->getNumComps() == 1) {
    std::cout << "Only 1 B-component, DivB matrix is zero!\n";
    MxCrsMatrix<Scalar>::fillComplete(bfield->getMap(), psifield->getMap());
  }
  else
    setMatrix();
}

#if 0
template<>
void MxYeeDeyMittraDivB<1>::setMatrix() {;}

template<>
void MxYeeDeyMittraDivB<2>::setMatrix() {

  size_t numBComps = bfield->getNumComps();

  if (numBComps == 1) {
    std::cout << "MxYeeDeyMittraDivB::setMatrix(): 2D TE simulation, divergence of B "
              << "is always zero. This matrix is zero.\n";
    //matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, psifield->getMap(), 0));
  }
  else {
    //matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, psifield->getMap(), 4));

    MxDimVector<double, 2> cellSize(grid->getCellSize());
    double idx = 1.0 / cellSize[0];
    double idy = 1.0 / cellSize[1];

    MxGridFieldIter<2> fieldIter(psifield);

    int row;
    double val;
    int cols[4];
    //double vals[4] = {1, -1, -1, 1};
    double vals[4] = {idx, -idx, idy, -idy};
    MxComplex cfactors[4];
    double factors[4];

    MxDimVector<int, 2> cell;

    // (Div B)_ij = (B_{x|i+1,j} - B_{x|ij}) / dx + (B_{y|ij+1} - B_{y|ij}) / dy
    for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
      cell = fieldIter.getCell();
      row = psifield->globCompIndx(0, cell);

      cell[0]++;
      cols[0] = bfield->globCompIndx(0, cell);
      cfactors[0] = bfield->getCompFactor(0, cell);
      cell[0]--;
      cols[1] = bfield->globCompIndx(0, cell);
      cfactors[1] = bfield->getCompFactor(0, cell);
      cell[1]++;
      cols[2] = bfield->globCompIndx(1, cell);
      cfactors[2] = bfield->getCompFactor(1, cell);
      cell[1]--;
      cols[3] = bfield->globCompIndx(1, cell);
      cfactors[3] = bfield->getCompFactor(1, cell);

      for (size_t i = 0; i < 4; i++) {
        if (this->imaginary)
          factors[i] = cfactors[i].imag();
        else
          factors[i] = cfactors[i].real();

        if (factors[i] != 0.0) {
          val = factors[i] * vals[i];
          //matrix->InsertGlobalValues(row, 1, &vals[i], &cols[i]);
          Epetra_CrsMatrix::InsertGlobalValues(row, 1, &val, &cols[i]);
        }
      }
    }
    //matrix->FillComplete(bfield->getMap(), psifield->getMap());
    Epetra_CrsMatrix::FillComplete(bfield->getMap(), psifield->getMap());
  }
}

template<>
void MxYeeDeyMittraDivB<3>::setMatrix() {
  //matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, psifield->getMap(), 6));

  MxDimVector<double, 3> cellSize(grid->getCellSize());
  double idx = 1 / cellSize[0];
  double idy = 1 / cellSize[1];
  double idz = 1 / cellSize[2];

  int row;
  double val;
  int cols[6];
  //double vals[4] = {1, -1, -1, 1};
  double vals[6] = {idx, -idx, idy, -idy, idz, -idz};
  MxComplex cfactors[6];
  double factors[6];

  MxDimVector<int, 3> cell;

  MxGridFieldIter<3> fieldIter(psifield);

  for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
    cell = fieldIter.getCell();

    row = psifield->globCompIndx(0, cell);

    for (int comp = 0; comp < 3; comp++) {
      cell[comp]++;
      cols[2 * comp] = bfield->globCompIndx(comp, cell);
      cfactors[2 * comp] = bfield->getCompFactor(comp, cell);
      cell[comp]--;
      cols[2 * comp + 1] = bfield->globCompIndx(comp, cell);
      cfactors[2 * comp + 1] = bfield->getCompFactor(comp, cell);
    }

    //std::cout << "Frac : Global index : Local index\n";
    //int proc, lid;
    for (size_t i = 0; i < 6; i++) {
      if (this->imaginary)
        factors[i] = cfactors[i].imag();
      else
        factors[i] = cfactors[i].real();

      if (factors[i] != 0.0) {
        //efield->getMap().RemoteIDList(1, &cols[i], &proc, &lid);
        //std::cout << "  " << fracs[i] << " : " << cols[i] << " : " << lid << "\n";
        //matrix->InsertGlobalValues(row, 1, &vals[i], &cols[i]);
        //Epetra_CrsMatrix::InsertGlobalValues(row, 1, &vals[i], &cols[i]);
        val = factors[i] * vals[i];
        this->InsertGlobalValues(row, 1, &val, &cols[i]);
      }
    }
    //std::cout << "\n";
  }
  //matrix->FillComplete(bfield->getMap(), psifield->getMap());
  //Epetra_CrsMatrix::FillComplete(bfield->getMap(), psifield->getMap());
  this->FillComplete(bfield->getMap(), psifield->getMap());

}
#endif

template<size_t DIM, typename Scalar>
void MxYeeDeyMittraDivB<DIM, Scalar>::setMatrix() {

  MxDimVector<double, DIM> cellSize(grid->getCellSize());
  
  Scalar vals[2*DIM];
  for (size_t i = 0; i < DIM; ++i) {
    MxUtil::convertScalar(1.0 / cellSize[i], vals[2*i]);
    MxUtil::convertScalar(-1.0 / cellSize[i], vals[2*i+1]);
  }

  MxIndex row;
  Scalar val, factor;
  MxIndex cols[2*DIM];
  MxComplex factors[2*DIM];

  MxDimVector<int, DIM> cell;

  MxGridFieldIter<DIM> fieldIter(psifield);

  for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
    cell = fieldIter.getCell();

    row = psifield->globCompIndx(0, cell);

    for (size_t comp = 0; comp < DIM; comp++) {
      cell[comp]++;
      cols[2 * comp] = bfield->globCompIndx(comp, cell);
      factors[2 * comp] = bfield->getCompFactor(comp, cell);
      cell[comp]--;
      cols[2 * comp + 1] = bfield->globCompIndx(comp, cell);
      factors[2 * comp + 1] = bfield->getCompFactor(comp, cell);
    }

    for (size_t i = 0; i < 2*DIM; i++) {
      if (factors[i] != ScalarTraits<MxComplex>::zero()) {
        MxUtil::convertScalar(factors[i], factor);
        val = factor * vals[i];
        MxCrsMatrix<Scalar>::insertRowValues(row, 1, &cols[i], &val);
      }
    }
  }
  MxCrsMatrix<Scalar>::fillComplete(bfield->getMap(), psifield->getMap());

}


template class MxYeeDeyMittraDivB<1, double>;
template class MxYeeDeyMittraDivB<2, double>;
template class MxYeeDeyMittraDivB<3, double>;
template class MxYeeDeyMittraDivB<1, MxComplex>;
template class MxYeeDeyMittraDivB<2, MxComplex>;
template class MxYeeDeyMittraDivB<3, MxComplex>;
