#include "MxYeeElecFieldBase.h"
#include "MxGridDomainIter.hpp"
#include "MxUtil.hpp"


template<size_t DIM>
void MxYeeElecFieldBase<DIM>::initYeeElecFieldBase(
const MxGrid<DIM> * aGrid, MxPolType polarization) {
//lBCs(PERIODIC), uBCs(PERIODIC),
//useLowerComps(true), useUpperComps(false), 
//gridRes(aGrid->getResolution()), phaseShifts(0.0) {
  pol = polarization;
  this->fieldName = "YeeElecFieldBase";
  this->grid = aGrid;
  this->mapSet = false;
  this->regionSet = false;

  if (DIM == 1) {
    this->numComps = 1;
    this->compCellCoords.resize(1); // automatically initialized to 0
    vecCompDirs.resize(1);

    if (pol == TE) {
      vecCompDirs[0][1] = 1;
      compDirs.push_back(1);
    }
    else if (pol == TM) {
      vecCompDirs[0][2] = 1;
      compDirs.push_back(2);
    }
    else {
      std::cout << "TEM or invalid polarization chosen for MxYeeElecFieldBase<1>. Using TE (E_y).\n";
      vecCompDirs[0][1] = 1;
    }
  }
  else if (DIM == 2) {
    double dx = this->grid->getCellSize()[0];
    double dy = this->grid->getCellSize()[1];
    if (pol == TE) {
      this->numComps = 2;
      this->compCellCoords.resize(2); // automatically initialized to 0

      // put in x-edge component
      this->compCellCoords[0][0] = 0.5 * dx;
      vecCompDirs.push_back(vec3<double>(1, 0, 0));
      compDirs.push_back(0);
      // put in y-edge component
      this->compCellCoords[1][1] = 0.5 * dy;
      vecCompDirs.push_back(vec3<double>(0, 1, 0));
      compDirs.push_back(1);
    }
    else if (pol == TM) {
      this->numComps = 1;
      this->compCellCoords.resize(1); // automatically initialized to 0
      // put in z-point component
      vecCompDirs.push_back(vec3<double>(0, 0, 1));
      compDirs.push_back(2);
    }
    else {
      std::cout << "MxYeeElecFieldBase: invalid polarization.";
      throw 1;
    }
  }
  else {
    this->numComps = 3;
    double dx = this->grid->getCellSize()[0];
    double dy = this->grid->getCellSize()[1];
    double dz = this->grid->getCellSize()[2];
    this->compCellCoords.resize(3);

    // put in x-edge component
    this->compCellCoords[0][0] = 0.5 * dx;
    vecCompDirs.push_back(vec3<double>(1, 0, 0));
    compDirs.push_back(0);
    // put in y-edge component
    this->compCellCoords[1][1] = 0.5 * dy;
    vecCompDirs.push_back(vec3<double>(0, 1, 0));
    compDirs.push_back(1);
    // put in z-edge component
    this->compCellCoords[2][2] = 0.5 * dz;
    vecCompDirs.push_back(vec3<double>(0, 0, 1));
    compDirs.push_back(2);
  }

  this->initBCs();
}

template<size_t DIM>
void MxYeeElecFieldBase<DIM>::setBCs(MxDimVector<MxBCType, DIM> lower,
MxDimVector<MxBCType, DIM> upper) {
  MxBCType lbc, ubc;
  MxDimVector<MxBCType, DIM> lbcs(lower), ubcs(upper);

  for (size_t comp = 0; comp < this->numComps; ++comp) {
    // bcs for each component default to what was passed in
    lbcs = lower;
    ubcs = upper;
    for (size_t i = 0; i < DIM; ++i) {
      lbc = lower[i];
      ubc = upper[i];

      // Treat the specialized boundary conditions such as PEC, PMC
      // lower bcs
      if (lbc == PEC) {
        // if the component direction is normal to the boundary
        if (compDirs[comp] == i)
          lbcs[i] = CONSTANT;
        else
          lbcs[i] = ZERO;
      }
      else if (lbc == PMC) {
        // if the component direction is normal to the boundary
        if (compDirs[comp] == i)
          lbcs[i] = ZERO;
        else
          lbcs[i] = CONSTANT;
      }

      // upper bcs
      if (ubc == PEC) {
        // if the component direction is normal to the boundary
        if (compDirs[comp] == i)
          ubcs[i] = CONSTANT;
        else
          ubcs[i] = ZERO;
      }
      else if (ubc == PMC) {
        // if the component direction is normal to the boundary
        if (compDirs[comp] == i)
          ubcs[i] = ZERO;
        else
          ubcs[i] = CONSTANT;
      }

    }
    this->setCompBCs(comp, lbcs, ubcs);
  }
}

// explicit specialization
template class MxYeeElecFieldBase<1>;
template class MxYeeElecFieldBase<2>;
template class MxYeeElecFieldBase<3>;
