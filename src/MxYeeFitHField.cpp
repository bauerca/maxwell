#include "MxYeeFitHField.h"
#include "MxGridDomainIter.hpp"
#include "MxUtil.hpp"


template<size_t DIM>
MxYeeFitHField<DIM>::MxYeeFitHField(const MxGrid<DIM> * aGrid,
MxPolType polarization) {
  this->initYeeMagFieldBase(aGrid, polarization);
  this->fieldName = "YeeFitHField";

  if (DIM == 1) {
    this->compPtopes.push_back(Teuchos::rcp(new MxPoint<DIM>()));
  }
  else if (DIM == 2) {
    double dx = this->grid->getCellSize()[0];
    double dy = this->grid->getCellSize()[1];
    if (polarization == TE) {
      // put in Hz-point component
      this->compPtopes.push_back(Teuchos::rcp(new MxPoint<DIM>()));
    }
    else if (polarization == TM) {
      // put in x-edge component
      this->compPtopes.push_back(Teuchos::rcp(new MxCartSeg<DIM>('x', dx)));
      // put in y-edge component
      this->compPtopes.push_back(Teuchos::rcp(new MxCartSeg<DIM>('y', dy)));
    }
    else {
      std::cout << "MxYeeFitHField: invalid polarization.";
      throw 1;
    }
  }
  else {
    double dx = this->grid->getCellSize()[0];
    double dy = this->grid->getCellSize()[1];
    double dz = this->grid->getCellSize()[2];

    // put in x-edge component
    this->compPtopes.push_back(Teuchos::rcp(new MxCartSeg<DIM>('x', dx)));
    // put in y-edge component
    this->compPtopes.push_back(Teuchos::rcp(new MxCartSeg<DIM>('y', dy)));
    // put in z-edge component
    this->compPtopes.push_back(Teuchos::rcp(new MxCartSeg<DIM>('z', dz)));
  }

}

#if 0
template<size_t DIM>
void MxYeeFitHField<DIM>::setBCs(Teuchos::ParameterList theBCList) {
  bcList = theBCList;

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

template<size_t DIM>
MxDimVector<int, DIM> MxYeeFitHField<DIM>::getInteriorComp(size_t comp, MxDimVector<int, DIM> cell) const {
  MxDimVector<int, DIM> newCell(cell);
  bool compDirMatchDir;
  for (size_t i = 0; i < DIM; ++i) {
    compDirMatchDir = (compDirs[comp] == i);
    if (cell[i] >= gridRes[i]) {
      switch (uBCs[i]) {
        case SYMM:
          newCell[i] += gridRes[i] - cell[i] - (compDirMatchDir ? 0 : 1);
          break;
        case ANTISYMM:
          newCell[i] += gridRes[i] - cell[i] - (compDirMatchDir ? 0 : 1);
          break;
        case PEC:
          newCell[i] = 2 * gridRes[i] - cell[i] - (compDirMatchDir ? 0 : 1);
          break;
        case PMC:
          newCell[i] = 2 * gridRes[i] - cell[i] - (compDirMatchDir ? 0 : 1);
          break;
        case PERIODIC:
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
        case PERIODIC:
          newCell[i] = gridRes[i] + cell[i];
      }
    }
  }

  return newCell;
}



template<size_t DIM>
size_t MxYeeFitHField<DIM>::globCompIndx(size_t comp, const MxDimVector<int, DIM> & cell) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  return comp + this->numComps * this->grid->cellToGlobalIndx(newCell);
}

template<size_t DIM>
double MxYeeFitHField<DIM>::calcCompFrac(size_t comp, MxDimVector<int, DIM> cell, const MxShape<DIM> & aShape) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  MxDimVector<double, DIM> p(this->grid->nodeCoord(newCell) + this->compCellCoords[comp]);
  return this->compPtopes[comp]->volumeFraction(aShape, p);
}

template<size_t DIM>
MxComplex MxYeeFitHField<DIM>::getCompFactor(size_t comp, MxDimVector<int, DIM> cell) const {
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


// replaces below, which has been put back in MxGridField base class
template<size_t DIM>
bool MxYeeFitHField<DIM>::useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const {
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
#endif

// explicit specialization
template class MxYeeFitHField<1>;
template class MxYeeFitHField<2>;
template class MxYeeFitHField<3>;
