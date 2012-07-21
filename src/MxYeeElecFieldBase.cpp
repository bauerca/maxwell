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

#if 0
template<size_t DIM>
void MxYeeElecFieldBase<DIM>::setBCs(Teuchos::ParameterList bcList) {

  lBCs = bcList.get("lower bcs", MxDimVector<MxBCType, DIM>(PERIODIC)); 
  uBCs = bcList.get("upper bcs", MxDimVector<MxBCType, DIM>(PERIODIC)); 
  phaseShifts = bcList.get("phase shifts", MxDimVector<double, DIM>(0.0));

  //std::cout << "upper bcs: " << uBCs;
  //std::cout << "lower bcs: " << lBCs;

  for (size_t i = 0; i < DIM; ++i) {
    //useUpperComps[i] = (uBCs[i] == SYMM) or (uBCs[i] == ANTISYMM);
    useUpperComps[i] = (uBCs[i] == PMC);
    useLowerComps[i] = (lBCs[i] != PEC); 
  }
}
#endif

// moved to super class (12/30/2011)
#if 0
template<size_t DIM>
MxDimVector<int, DIM> MxYeeElecFieldBase<DIM>::getInteriorComp(size_t comp, MxDimVector<int, DIM> cell) const {
  MxDimVector<int, DIM> newCell(cell);
  bool compDirMatchDir;
  for (size_t i = 0; i < DIM; ++i) {
    compDirMatchDir = (compDirs[comp] == i);
    if (cell[i] >= gridRes[i]) {
      switch (uBCs[i]) {
        case SYMM:
          newCell[i] += gridRes[i] - cell[i] - (compDirMatchDir ? 1 : 0);
          break;
        case ANTISYMM:
          newCell[i] += gridRes[i] - cell[i] - (compDirMatchDir ? 1 : 0);
          break;
        case PEC:
          newCell[i] = 2 * gridRes[i] - cell[i] - (compDirMatchDir ? 1 : 0);
          break;
        case PMC:
          newCell[i] = 2 * gridRes[i] - cell[i] - (compDirMatchDir ? 1 : 0);
          break;
        case PERIODIC:
          newCell[i] = cell[i] - gridRes[i];
      }
    }
    else if (cell[i] < 0) {
      switch (lBCs[i]) {
        case PEC:
          newCell[i] = -cell[i] - (compDirMatchDir ? 1 : 0);
          break;
        case PMC:
          newCell[i] = -cell[i] - (compDirMatchDir ? 1 : 0);
          break;
        case PERIODIC:
          newCell[i] = gridRes[i] + cell[i];
      }
    }
  }

  return newCell;
}



template<size_t DIM>
size_t MxYeeElecFieldBase<DIM>::globCompIndx(size_t comp, const MxDimVector<int, DIM> & cell) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  return comp + this->numComps * this->grid->cellToGlobalIndx(newCell);
}

template<size_t DIM>
double MxYeeElecFieldBase<DIM>::calcCompFrac(size_t comp, MxDimVector<int, DIM> cell, const MxShape<DIM> & aShape) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  MxDimVector<double, DIM> p(this->grid->nodeCoord(newCell) + this->compCellCoords[comp]);
  return this->compPtopes[comp]->volumeFraction(aShape, p);
}

template<size_t DIM>
MxComplex MxYeeElecFieldBase<DIM>::getCompFactor(size_t comp, MxDimVector<int, DIM> cell) const {
  if (this->regionSet)
    if (this->getCompFrac(comp, cell, this->regionName) == 0.0)
      return 0.0;

  // compInMap interiorizes the component, then checks if it exists
  // if interiorized component does not exist in map, then its factor must be zero
  // cab: This is so sloooww! e.g. when forming curl matrix, which calls getCompFactor
  //if (!this->compInMap(comp, cell))
  //  return 0.0;

  MxComplex res = 1.0;
  bool compDirMatchDir;
  for (size_t i = 0; i < DIM; ++i) {
    compDirMatchDir = (compDirs[comp] == i);
    if (cell[i] == gridRes[i]) {
      switch (uBCs[i]) {
        case ANTISYMM:
          res *= (compDirMatchDir ? -1.0 : 1.0);
          break;
        case PEC:
          res *= (compDirMatchDir ? 1.0 : 0.0);
          break;
        case PMC:
          res *= (compDirMatchDir ? -1.0 : 1.0);
          break;
        default:
          res *= 1.0;
      }
    }
    else if (cell[i] == 0) {
      switch (lBCs[i]) {
        case ANTISYMM:
          res *= (compDirMatchDir ? -1.0 : 1.0);
          break;
        case PEC:
          res *= (compDirMatchDir ? 1.0 : 0.0);
          break;
        default:
          res *= 1.0;
      }
    }
    else if (cell[i] > gridRes[i]) {
      switch (uBCs[i]) {
        case ANTISYMM:
          res *= -1.0;
          break;
        case PEC:
          res *= (compDirMatchDir ? 1.0 : -1.0);
          break;
        case PMC:
          res *= (compDirMatchDir ? -1.0 : 1.0);
          break;
        default:
          res *= 1.0;
      }
    }
    else if (cell[i] < 0) {
      switch (lBCs[i]) {
        case ANTISYMM:
          res *= -1.0;
          break;
        case PEC:
          res *= (compDirMatchDir ? 1.0 : -1.0);
          break;
        case PMC:
          res *= (compDirMatchDir ? -1.0 : 1.0);
          break;
        default:
          res *= 1.0;
      }
    }
  }

  // sign is set, now set mag with phase shift
  MxDimVector<double, DIM> coord(this->getCompCoord(comp, cell));
  MxDimVector<double, DIM> ub(this->grid->getUpperBounds());
  MxDimVector<double, DIM> lb(this->grid->getLowerBounds());
  MxComplex I(0.0, 1.0);
  for (size_t i = 0; i < DIM; ++i) {
    if (lBCs[i] == PERIODIC && uBCs[i] == PERIODIC) {
      if (coord[i] >= ub[i])
        res *= exp(I * phaseShifts[i]);
      else if (coord[i] < lb[i])
        res *= exp(-I * phaseShifts[i]);
    }
  }

  return res;
}


// replaces below, which has been put back in MxGridField base class
template<size_t DIM>
bool MxYeeElecFieldBase<DIM>::useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const {
  if (this->regionSet and this->getCompFrac(comp, cell, this->regionName) == 0.0) return false;

  for (size_t i = 0; i < DIM; ++i) {
    if (cell[i] == gridRes[i] and (compDirs[comp] == i or !useUpperComps[i])) {
      return false;
    }
    else if (cell[i] == 0 and compDirs[comp] != i and !useLowerComps[i]) {
      return false;
    }
  }
  return true;
}
#endif // moved to base class


// explicit specialization
template class MxYeeElecFieldBase<1>;
template class MxYeeElecFieldBase<2>;
template class MxYeeElecFieldBase<3>;
