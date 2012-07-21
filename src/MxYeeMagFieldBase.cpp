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

// moved to super class (12/30/2011)
#if 0
template<size_t DIM>
void MxYeeMagField<DIM>::setBCs(Teuchos::ParameterList theBCList) {
  bcList = theBCList;

  lBCs = bcList.get("lower bcs", MxDimVector<MxBCType, DIM>(PERIODIC)); 
  uBCs = bcList.get("upper bcs", MxDimVector<MxBCType, DIM>(PERIODIC)); 
  phaseShifts = bcList.get("phase shifts", MxDimVector<double, DIM>(0));

  for (size_t i = 0; i < DIM; ++i) {
    //useUpperComps[i] = (uBCs[i] == SYMM) or (uBCs[i] == ANTISYMM);
    useUpperComps[i] = (uBCs[i] == PMC);
    useLowerComps[i] = (lBCs[i] != PEC);
  }
}


template<size_t DIM>
MxDimVector<int, DIM> MxYeeMagField<DIM>::getInteriorComp(size_t comp, MxDimVector<int, DIM> cell) const {
  MxDimVector<int, DIM> newCell(cell);
  bool compDirMatchDir;
  for (size_t i = 0; i < DIM; ++i) {
    compDirMatchDir = (compDirs[comp] == i);
    if (cell[i] >= gridRes[i]) {
      switch (uBCs[i]) {
        case SYMM:
          newCell[i] = 2 * gridRes[i] - cell[i] - (compDirMatchDir ? 0 : 1);
          break;
        case ANTISYMM:
          newCell[i] = 2 * gridRes[i] - cell[i] - (compDirMatchDir ? 0 : 1);
          break;
        case PEC:
          newCell[i] = 2 * gridRes[i] - cell[i] - (compDirMatchDir ? 0 : 1);
          break;
        case PMC:
          newCell[i] = 2 * gridRes[i] - cell[i] - (compDirMatchDir ? 0 : 1);
          break;
        default:
          newCell[i] = cell[i] - gridRes[i];
      }
    }
    else if (cell[i] < 0) {
      switch (lBCs[i]) {
        case PEC:
          newCell[i] = -cell[i] - (compDirMatchDir ? 0 : 1);
          break;
        case PMC:
          newCell[i] = -cell[i] - (compDirMatchDir ? 0 : 1);
          break;
        default:
          newCell[i] = gridRes[i] + cell[i];
      }
    }
  }

  return newCell;
}



template<size_t DIM>
size_t MxYeeMagField<DIM>::globCompIndx(size_t comp, const MxDimVector<int, DIM> & cell) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  return comp + this->numComps * this->grid->cellToGlobalIndx(newCell);
}

template<size_t DIM>
double MxYeeMagField<DIM>::calcCompFrac(size_t comp, MxDimVector<int, DIM> cell, const MxShape<DIM> & aShape) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  MxDimVector<double, DIM> p(this->grid->nodeCoord(newCell) + this->compCellCoords[comp]);
  return this->compPtopes[comp]->volumeFraction(aShape, p);
}


template<size_t DIM>
MxComplex MxYeeMagField<DIM>::getCompFactor(size_t comp,
MxDimVector<int, DIM> cell) const {
  if (this->regionSet)
    if (this->getCompFrac(comp, cell, this->regionName) == 0.0)
      return 0.0;

  // compInMap interiorizes the component, then checks if it exists
  //if (!this->compInMap(comp, cell))
  //  return 0.0;

  MxComplex res = 1.0;
  bool compDirMatchDir;
  for (size_t i = 0; i < DIM; ++i) {
    compDirMatchDir = (compDirs[comp] == i);
    if (cell[i] == gridRes[i]) {
      switch (uBCs[i]) {
        case PEC:
          res *= (compDirMatchDir ? 0.0 : 1.0);
          break;
        case PMC:
          res *= (compDirMatchDir ? 1.0 : -1.0);
          break;
        default:
          res *= 1.0;
      }
    }
    else if (cell[i] == 0) {
      switch (lBCs[i]) {
        case PEC:
          res *= (compDirMatchDir ? 0.0 : 1.0);
          break;
        default:
          res *= 1.0;
      }
    }
    else if (cell[i] > gridRes[i]) {
      switch (uBCs[i]) {
        case PEC:
          res *= (compDirMatchDir ? -1.0 : 1.0);
          break;
        case PMC:
          res *= (compDirMatchDir ? 1.0 : -1.0);
          break;
        default:
          res *= 1.0;
      }
    }
    else if (cell[i] < 0) {
      switch (lBCs[i]) {
        case PEC:
          res *= (compDirMatchDir ? -1.0 : 1.0);
          break;
        case PMC:
          res *= (compDirMatchDir ? 1.0 : -1.0);
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


// replaces below. This function is called by MxGridField::setEpetraMap()
template<size_t DIM>
bool MxYeeMagField<DIM>::useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const {
  if (this->regionSet and this->getCompFrac(comp, cell, this->regionName) == 0.0) return false;

  for (size_t i = 0; i < DIM; ++i) {
    if (cell[i] == gridRes[i] and (compDirs[comp] != i or !useUpperComps[i])) {
      return false;
    }
    else if (cell[i] == 0 and compDirs[comp] == i and !useLowerComps[i]) {
      return false;
    }
  }
  return true;
}


// overloaded getInterpolator
template<size_t DIM>
Epetra_CrsMatrix MxYeeMagField<DIM>::getInterpolator(const MxGridField<DIM> & targetField,
bool topeIntegrated) const {
  Epetra_CrsMatrix res(MxGridField<DIM>::getInterpolator(targetField));

  if (topeIntegrated) {
    if (this->hasShapeRep("pec") and targetField.hasShapeRep("pec")) {
      Epetra_Vector srcFracs(this->getAllCompFracs("pec"));
      srcFracs.Reciprocal(srcFracs);
      res.RightScale(srcFracs);
      res.LeftScale(targetField.getAllCompFracs("pec"));
    }
  }

  return res;
}

template<size_t DIM>
Epetra_CrsMatrix MxYeeMagField<DIM>::getInterpolator(const MxGridField<DIM> & targetField) const {
  return MxGridField<DIM>::getInterpolator(targetField);
}

#endif


// explicit specialization
template class MxYeeMagFieldBase<1>;
template class MxYeeMagFieldBase<2>;
template class MxYeeMagFieldBase<3>;
