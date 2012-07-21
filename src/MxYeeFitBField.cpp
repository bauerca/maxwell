#include "MxYeeFitBField.h"
#include "MxPoint.hpp"
#include "MxCartSeg.hpp"
#include "MxCartRect.hpp"
#include "MxGridDomainIter.hpp"
#include "MxUtil.hpp"


template<>
MxYeeFitBField<1>::MxYeeFitBField(const MxGrid<1> * aGrid,
MxPolType polarization) : mDMFrac(0) {
  this->initYeeMagFieldBase(aGrid, polarization);
  this->fieldName = "YeeFitBField";

  this->compPtopes.push_back(Teuchos::rcp(new MxPoint<1>()));
}

template<>
MxYeeFitBField<2>::MxYeeFitBField(const MxGrid<2> * aGrid,
MxPolType polarization) : mDMFrac(0) {
  this->initYeeMagFieldBase(aGrid, polarization);
  this->fieldName = "YeeFitBField";

  double dx = this->grid->getCellSize()[0];
  double dy = this->grid->getCellSize()[1];
  if (polarization == TM) {
    // put in x-face component
    this->compPtopes.push_back(Teuchos::rcp(new MxCartSeg<2>('y', dy)));
    // put in y-face component
    this->compPtopes.push_back(Teuchos::rcp(new MxCartSeg<2>('x', dx)));
  }
  else if (polarization == TE) {
    // put in z-face component
    this->compPtopes.push_back(Teuchos::rcp(new MxCartRect<2>('z', dx, dy)));
  }
  else {
    std::cout << "MxYeeFitBField: invalid polarization.";
    throw 1;
  }
}


template<>
MxYeeFitBField<3>::MxYeeFitBField(const MxGrid<3> * aGrid,
MxPolType polarization) : mDMFrac(0) {
  this->initYeeMagFieldBase(aGrid, polarization);
  this->fieldName = "YeeFitBField";

  double dx = this->grid->getCellSize()[0];
  double dy = this->grid->getCellSize()[1];
  double dz = this->grid->getCellSize()[2];

  // put in x-face component
  this->compPtopes.push_back(Teuchos::rcp(new MxCartRect<3>('x', dy, dz)));
  // put in y-face component
  this->compPtopes.push_back(Teuchos::rcp(new MxCartRect<3>('y', dz, dx)));
  // put in z-face component
  this->compPtopes.push_back(Teuchos::rcp(new MxCartRect<3>('z', dx, dy)));
} 

// for Dey-Mittra
template<size_t DIM>
bool MxYeeFitBField<DIM>::useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const {
  bool res = MxYeeMagFieldBase<DIM>::useCompInMap(comp, cell);

  if (this->regionSet) {
    if (this->getCompFrac(comp, cell, this->regionName) < mDMFrac) {
      res = false;
    }
  }

  return res;
}

template<size_t DIM>
MxComplex MxYeeFitBField<DIM>::getCompFactor(size_t comp,
MxDimVector<int, DIM> cell) const {
  if (this->regionSet)
    if (this->getCompFrac(comp, cell, this->regionName) < mDMFrac)
      return 0.0;

  return MxYeeMagFieldBase<DIM>::getCompFactor(comp, cell);
}

#if 0
template<size_t DIM>
void MxYeeFitBField<DIM>::setBCs(Teuchos::ParameterList theBCList) {
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
MxDimVector<int, DIM> MxYeeFitBField<DIM>::getInteriorComp(size_t comp, MxDimVector<int, DIM> cell) const {
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
size_t MxYeeFitBField<DIM>::globCompIndx(size_t comp, const MxDimVector<int, DIM> & cell) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  return comp + this->numComps * this->grid->cellToGlobalIndx(newCell);
}

template<size_t DIM>
double MxYeeFitBField<DIM>::calcCompFrac(size_t comp, MxDimVector<int, DIM> cell, const MxShape<DIM> & aShape) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  MxDimVector<double, DIM> p(this->grid->nodeCoord(newCell) + this->compCellCoords[comp]);
  return this->compPtopes[comp]->volumeFraction(aShape, p);
}


template<size_t DIM>
MxComplex MxYeeFitBField<DIM>::getCompFactor(size_t comp,
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
bool MxYeeFitBField<DIM>::useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const {
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
Epetra_CrsMatrix MxYeeFitBField<DIM>::getInterpolator(const MxGridField<DIM> & targetField,
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
Epetra_CrsMatrix MxYeeFitBField<DIM>::getInterpolator(const MxGridField<DIM> & targetField) const {
  return MxGridField<DIM>::getInterpolator(targetField);
}

#endif


// explicit specialization
template class MxYeeFitBField<1>;
template class MxYeeFitBField<2>;
template class MxYeeFitBField<3>;
