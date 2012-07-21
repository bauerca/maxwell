#include "MxYeePsiField.h"
#include "MxCartSeg.hpp"
#include "MxCartRect.hpp"
#include "MxCartBox.hpp"
#include "MxGridDomainIter.hpp"
#include "MxUtil.hpp"


template<>
MxYeePsiField<1>::MxYeePsiField(const MxGrid<1> * aGrid, const MxYeeFitBField<1> * bfield) :
mBfield(bfield) {
  this->grid = aGrid;
  this->mapSet = false;
  this->regionSet = false;
  this->fieldName = "YeeFitPsiField";

  double dx = this->grid->getCellSize()[0];

  this->numComps = 1;
  this->compCellCoords.resize(1); // automatically initialized to 0

  this->compPtopes.push_back(Teuchos::rcp(new MxCartSeg<1>('x', dx)));
  this->compCellCoords[0][0] = 0.5 * dx;

  this->initBCs();
}

template<>
MxYeePsiField<2>::MxYeePsiField(const MxGrid<2> * aGrid, const MxYeeFitBField<2> * bfield) :
mBfield(bfield) {
  this->grid = aGrid;
  this->mapSet = false;
  this->regionSet = false;
  this->fieldName = "YeeFitPsiField";

  double dx = this->grid->getCellSize()[0];
  double dy = this->grid->getCellSize()[1];

  this->numComps = 1;
  this->compCellCoords.resize(1); // automatically initialized to 0

  this->compPtopes.push_back(Teuchos::rcp(new MxCartRect<2>('z', dx, dy)));
  this->compCellCoords[0][0] = 0.5 * dx;
  this->compCellCoords[0][1] = 0.5 * dy;

  this->initBCs();
}


template<>
MxYeePsiField<3>::MxYeePsiField(const MxGrid<3> * aGrid, const MxYeeFitBField<3> * bfield) : 
mBfield(bfield) {
  this->grid = aGrid;
  this->mapSet = false;
  this->regionSet = false;
  this->fieldName = "YeeFitPsiField";

  this->numComps = 1;
  double dx = this->grid->getCellSize()[0];
  double dy = this->grid->getCellSize()[1];
  double dz = this->grid->getCellSize()[2];
  this->compCellCoords.resize(1);

  this->compPtopes.push_back(Teuchos::rcp(new MxCartBox(dx, dy, dz)));
  this->compCellCoords[0][0] = 0.5 * dx;
  this->compCellCoords[0][1] = 0.5 * dy;
  this->compCellCoords[0][2] = 0.5 * dz;

  this->initBCs();
} 

template<size_t DIM>
void MxYeePsiField<DIM>::setBCs(MxDimVector<MxBCType, DIM> lower,
MxDimVector<MxBCType, DIM> upper) {
  MxBCType lbc, ubc;
  MxDimVector<MxBCType, DIM> lbcs(lower), ubcs(upper);

  for (size_t i = 0; i < DIM; ++i) {
    lbc = lower[i];
    ubc = upper[i];

    // Treat the specialized boundary conditions such as PEC, PMC
    // lower bcs
    if (lbc == PEC)
      lbcs[i] = CONSTANT; // so that normal gradient of Psi = 0
    else if (lbc == PMC)
      lbcs[i] = ZERO;

    // upper bcs
    if (ubc == PEC)
      ubcs[i] = CONSTANT; // so that normal gradient of Psi = 0
    else if (ubc == PMC)
      ubcs[i] = ZERO;
  }
  this->setCompBCs(0, lbcs, ubcs);
}


template<size_t DIM>
bool MxYeePsiField<DIM>::useCompInMap(size_t comp,
MxDimVector<int, DIM> cell) const {

  bool res = MxGridField<DIM>::useCompInMap(comp, cell);
  if (res == false) return res;

  for (int i = 0; i < DIM; ++i) {
    if (mBfield->useCompInMap(i, cell)) return true;
    cell[i]++;
    if (mBfield->useCompInMap(i, cell)) return true;
    cell[i]--;
  }

  return false;
}

template<size_t DIM>
MxComplex MxYeePsiField<DIM>::getCompFactor(size_t comp,
MxDimVector<int, DIM> cell) const {

  MxComplex res = MxGridField<DIM>::getCompFactor(comp, cell);

  if (!this->useCompInMap(comp, cell)) return 0.0;
  else return res;
}

#if 0
template<size_t DIM>
void MxYeePsiField<DIM>::setBCs(Teuchos::ParameterList theBCList) {
  bcList = theBCList;

  lBCs = bcList.get("lower bcs", MxDimVector<MxBCType, DIM>(PERIODIC)); 
  uBCs = bcList.get("upper bcs", MxDimVector<MxBCType, DIM>(PERIODIC)); 
  phaseShifts = bcList.get("phase shifts", MxDimVector<double, DIM>(0.0));
}

template<size_t DIM>
MxDimVector<int, DIM> MxYeePsiField<DIM>::getInteriorComp(size_t comp, MxDimVector<int, DIM> cell) const {
  MxDimVector<int, DIM> newCell(cell);
  for (size_t i = 0; i < DIM; ++i) {
    if (cell[i] >= gridRes[i]) {
      switch (uBCs[i]) {
        case PEC:
          newCell[i] = 2 * gridRes[i] - cell[i] - 1;
          break;
        case PMC:
          newCell[i] = 2 * gridRes[i] - cell[i] - 1;
          break;
        default:
          newCell[i] = cell[i] - gridRes[i];
      }
    }
    else if (cell[i] < 0) {
      switch (lBCs[i]) {
        case PEC:
          newCell[i] = -cell[i] - 1;
          break;
        case PMC:
          newCell[i] = -cell[i] - 1;
          break;
        default:
          newCell[i] = gridRes[i] + cell[i];
      }
    }
  }

  return newCell;
}

template<size_t DIM>
size_t MxYeePsiField<DIM>::globCompIndx(size_t comp, const MxDimVector<int, DIM> & cell) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  return comp + this->numComps * this->grid->cellToGlobalIndx(newCell);
}

template<size_t DIM>
double MxYeePsiField<DIM>::calcCompFrac(size_t comp, MxDimVector<int, DIM> cell, const MxShape<DIM> & aShape) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  MxDimVector<double, DIM> p(this->grid->nodeCoord(newCell) + this->compCellCoords[comp]);
  return this->compPtopes[comp]->volumeFraction(aShape, p);
}


template<size_t DIM>
MxComplex MxYeePsiField<DIM>::getCompFactor(size_t comp,
MxDimVector<int, DIM> cell) const {
  if (this->regionSet)
    if (this->getCompFrac(comp, cell, this->regionName) == 0.0)
      return 0.0;
  
  // compInMap interiorizes the component, then checks if it exists
  //if (!this->compInMap(comp, cell))
  //  return 0.0;

  MxComplex res = 1.0;
  for (size_t i = 0; i < DIM; ++i) {
    if (cell[i] >= gridRes[i]) {
      switch (uBCs[i]) {
        case ANTISYMM:
          res *= -1.0;
          break;
        case PMC:
          res *= -1.0;
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
        case PMC:
          res *= -1.0;
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

template<size_t DIM>
bool MxYeePsiField<DIM>::useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const {
  if (this->regionSet and this->getCompFrac(comp, cell, this->regionName) == 0.0) return false;

  for (size_t i = 0; i < DIM; ++i)
    if (cell[i] >= gridRes[i] or cell[i] < 0) return false;
  return true;
}
          
#endif

// explicit specialization
template class MxYeePsiField<1>;
template class MxYeePsiField<2>;
template class MxYeePsiField<3>;
