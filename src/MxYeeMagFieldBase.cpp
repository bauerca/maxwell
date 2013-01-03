#include "MxYeeMagFieldBase.h"
#include "MxPoint.hpp"
#include "MxCartSeg.hpp"
#include "MxCartRect.hpp"
#include "MxGridDomainIter.hpp"
#include "MxUtil.hpp"


template<>
void MxYeeMagFieldBase<1>::initYeeMagFieldBase(const MxGrid<1> * aGrid,
MxPolType polarization) {
  pol = polarization;
  this->fieldName = "YeeMagFieldBase";
  this->grid = aGrid;
  this->mapSet = false;
  this->regionSet = false;

  this->numComps = 1;
  this->compCellCoords.resize(1); // automatically initialized to 0
  vecCompDirs.resize(1);
  compDirs.resize(1);

  if (pol == TE) {
    vecCompDirs[0][1] = 1;
    compDirs[0] = 1;
  }
  else if (pol == TM) {
    vecCompDirs[0][2] = 1;
    compDirs[0] = 2;
  }
  else {
    std::cout << "invalid polarization chosen for MxYeeMagFieldBase<1>. Using TE (E_y).\n";
    vecCompDirs[0][1] = 1;
    compDirs[0] = 1;
  }

  this->initBCs();

}

template<>
void MxYeeMagFieldBase<2>::initYeeMagFieldBase(const MxGrid<2> * aGrid,
MxPolType polarization) {
  pol = polarization;
  this->fieldName = "YeeMagFieldBase";
  this->grid = aGrid;
  this->mapSet = false;
  this->regionSet = false;

  double dx = this->grid->getCellSize()[0];
  double dy = this->grid->getCellSize()[1];
  if (pol == TM) {
    this->numComps = 2;
    this->compCellCoords.resize(2); // automatically initialized to 0

    // put in x-face component location
    this->compCellCoords[0][1] = 0.5 * dy;
    vecCompDirs.push_back(vec3<double>(1, 0, 0));
    compDirs.push_back(0);
    // put in y-face component location
    this->compCellCoords[1][0] = 0.5 * dx;
    vecCompDirs.push_back(vec3<double>(0, 1, 0));
    compDirs.push_back(1);
  }
  else if (pol == TE) {
    this->numComps = 1;
    this->compCellCoords.resize(1); // automatically initialized to 0
    // put in z-face component location
    this->compCellCoords[0][0] = 0.5 * dx;
    this->compCellCoords[0][1] = 0.5 * dy;
    vecCompDirs.push_back(vec3<double>(0, 0, 1));
    compDirs.push_back(2);
  }
  else {
    std::cout << "MxYeeMagField: invalid polarization.";
    throw 1;
  }

  this->initBCs();
}
  


template<>
void MxYeeMagFieldBase<3>::initYeeMagFieldBase(const MxGrid<3> * aGrid,
MxPolType polarization) {
  pol = polarization;
  this->fieldName = "YeeMagField";
  this->grid = aGrid;
  this->mapSet = false;
  this->regionSet = false;

  this->numComps = 3;
  double dx = this->grid->getCellSize()[0];
  double dy = this->grid->getCellSize()[1];
  double dz = this->grid->getCellSize()[2];
  this->compCellCoords.resize(3);

  // put in x-face component location
  this->compCellCoords[0][1] = 0.5 * dy;
  this->compCellCoords[0][2] = 0.5 * dz;
  vecCompDirs.push_back(vec3<double>(1, 0, 0));
  compDirs.push_back(0);
  // put in y-face component location
  this->compCellCoords[1][2] = 0.5 * dz;
  this->compCellCoords[1][0] = 0.5 * dx;
  vecCompDirs.push_back(vec3<double>(0, 1, 0));
  compDirs.push_back(1);
  // put in z-face component location
  this->compCellCoords[2][0] = 0.5 * dx;
  this->compCellCoords[2][1] = 0.5 * dy;
  vecCompDirs.push_back(vec3<double>(0, 0, 1));
  compDirs.push_back(2);

  this->initBCs();
} 

template<size_t DIM>
void MxYeeMagFieldBase<DIM>::setBCs(MxDimVector<MxBCType, DIM> lower,
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
          lbcs[i] = ZERO;
        else
          lbcs[i] = CONSTANT;
      }
      else if (lbc == PMC) {
        // if the component direction is normal to the boundary
        if (compDirs[comp] == i)
          lbcs[i] = CONSTANT;
        else
          lbcs[i] = ZERO;
      }

      // upper bcs
      if (ubc == PEC) {
        // if the component direction is normal to the boundary
        if (compDirs[comp] == i)
          ubcs[i] = ZERO;
        else
          ubcs[i] = CONSTANT;
      }
      else if (ubc == PMC) {
        // if the component direction is normal to the boundary
        if (compDirs[comp] == i)
          ubcs[i] = CONSTANT;
        else
          ubcs[i] = ZERO;
      }

    }
    this->setCompBCs(comp, lbcs, ubcs);
  }
}



// explicit specialization
template class MxYeeMagFieldBase<1>;
template class MxYeeMagFieldBase<2>;
template class MxYeeMagFieldBase<3>;
