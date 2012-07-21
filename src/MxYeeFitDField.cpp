#include "MxYeeFitDField.h"
#include "MxPoint.hpp"
#include "MxCartSeg.hpp"
#include "MxCartRect.hpp"
#include "MxGridDomainIter.hpp"
#include "MxUtil.hpp"


template<>
MxYeeFitDField<1>::MxYeeFitDField(const MxGrid<1> * aGrid,
MxPolType polarization) {
  this->initYeeElecFieldBase(aGrid, polarization);
  this->fieldName = "yeeFitDField";

  this->compPtopes.push_back(Teuchos::rcp(new MxPoint<1>()));
}

template<>
MxYeeFitDField<2>::MxYeeFitDField(const MxGrid<2> * aGrid,
MxPolType polarization) {
  this->initYeeElecFieldBase(aGrid, polarization);
  this->fieldName = "yeeFitDField";

  double dx = this->grid->getCellSize()[0];
  double dy = this->grid->getCellSize()[1];
  if (polarization == TE) {
    // put in x-dual-face component
    this->compPtopes.push_back(Teuchos::rcp(new MxCartSeg<2>('y', dy)));
    // put in y-dual-face component
    this->compPtopes.push_back(Teuchos::rcp(new MxCartSeg<2>('x', dx)));
  }
  else if (polarization == TM) {
    // put in z-point component
    this->compPtopes.push_back(Teuchos::rcp(new MxCartRect<2>('z', dx, dy)));
  }
  else {
    std::cout << "MxYeeFitDField: invalid polarization.";
    throw 1;
  }
}


template<>
MxYeeFitDField<3>::MxYeeFitDField(const MxGrid<3> * aGrid,
MxPolType polarization) {
  this->initYeeElecFieldBase(aGrid, polarization);
  this->fieldName = "yeeFitDField";

  double dx = this->grid->getCellSize()[0];
  double dy = this->grid->getCellSize()[1];
  double dz = this->grid->getCellSize()[2];

  // put in x-dual-face component
  this->compPtopes.push_back(Teuchos::rcp(new MxCartRect<3>('x', dy, dz)));
  // put in y-dual-face component
  this->compPtopes.push_back(Teuchos::rcp(new MxCartRect<3>('y', dz, dx)));
  // put in z-dual-face component
  this->compPtopes.push_back(Teuchos::rcp(new MxCartRect<3>('z', dx, dy)));
} 

#if 0
template<size_t DIM>
void MxYeeFitDField<DIM>::setBCs(Teuchos::ParameterList theBCList) {
  bcList = theBCList;

  lBCs = bcList.get("lower bcs", MxDimVector<MxBCType, DIM>(PERIODIC)); 
  uBCs = bcList.get("upper bcs", MxDimVector<MxBCType, DIM>(PERIODIC)); 
  phaseShifts = bcList.get("phase shifts", MxDimVector<double, DIM>(0.0));

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
MxDimVector<int, DIM> MxYeeFitDField<DIM>::getInteriorComp(size_t comp, MxDimVector<int, DIM> cell) const {
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
size_t MxYeeFitDField<DIM>::globCompIndx(size_t comp, const MxDimVector<int, DIM> & cell) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  return comp + this->numComps * this->grid->cellToGlobalIndx(newCell);
}

template<size_t DIM>
double MxYeeFitDField<DIM>::calcCompFrac(size_t comp, MxDimVector<int, DIM> cell, const MxShape<DIM> & aShape) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  MxDimVector<double, DIM> p(this->grid->nodeCoord(newCell) + this->compCellCoords[comp]);
  return this->compPtopes[comp]->volumeFraction(aShape, p);
}

template<size_t DIM>
MxComplex MxYeeFitDField<DIM>::getCompFactor(size_t comp,
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
bool MxYeeFitDField<DIM>::useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const {
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






// overloaded getInterpolator
template<size_t DIM>
Epetra_CrsMatrix MxYeeFitDField<DIM>::getInterpolator(const MxGridField<DIM> & targetField,
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
Epetra_CrsMatrix MxYeeFitDField<DIM>::getInterpolator(const MxGridField<DIM> & targetField) const {
  return MxGridField<DIM>::getInterpolator(targetField);
}
#endif


// explicit specialization
template class MxYeeFitDField<1>;
template class MxYeeFitDField<2>;
template class MxYeeFitDField<3>;
