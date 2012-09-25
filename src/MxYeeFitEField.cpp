#include "MxYeeFitEField.h"
#include "MxGridDomainIter.hpp"
#include "MxUtil.hpp"


template<size_t DIM>
MxYeeFitEField<DIM>::MxYeeFitEField(const MxGrid<DIM> * aGrid,
const MxYeeFitBField<DIM> * bfield,
MxPolType polarization) : mBfield(bfield) {
//lBCs(PERIODIC), uBCs(PERIODIC),
//useLowerComps(true), useUpperComps(false), 
//gridRes(aGrid->getResolution()), phaseShifts(0.0) {
  this->initYeeElecFieldBase(aGrid, polarization);
  this->fieldName = "YeeFitEField";

  if (DIM == 1) {
    this->compPtopes.push_back(Teuchos::rcp(new MxPoint<DIM>()));
  }
  else if (DIM == 2) {
    double dx = this->grid->getCellSize()[0];
    double dy = this->grid->getCellSize()[1];
    if (polarization == TE) {
      // put in x-edge component
      this->compPtopes.push_back(Teuchos::rcp(new MxCartSeg<DIM>('x', dx)));
      // put in y-edge component
      this->compPtopes.push_back(Teuchos::rcp(new MxCartSeg<DIM>('y', dy)));
    }
    else if (polarization == TM) {
      // put in z-point component
      this->compPtopes.push_back(Teuchos::rcp(new MxPoint<DIM>()));
    }
    else {
      std::cout << "MxYeeFitEField: invalid polarization.";
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

template<size_t DIM>
bool MxYeeFitEField<DIM>::useCompInMap(size_t comp,
MxDimVector<int, DIM> cell) const {

  bool res = MxYeeElecFieldBase<DIM>::useCompInMap(comp, cell);
  if (res == false) return res;

  if (DIM == 3 or (DIM == 2 and this->pol == TM)) {
    size_t comp2, comp3;
    if (DIM == 3) {
      comp2 = (comp + 1) % 3;
      comp3 = (comp + 2) % 3;
    } else {
      comp2 = 0;
      comp3 = 1;
    }

    if (!mBfield->useCompInMap(comp2, cell)) return false;
    cell[comp3]--;
    if (!mBfield->useCompInMap(comp2, cell)) return false;
    cell[comp3]++;

    if (!mBfield->useCompInMap(comp3, cell)) return false;
    cell[comp2]--;
    if (!mBfield->useCompInMap(comp3, cell)) return false;
    cell[comp2]++;
  } else if (DIM == 2) {
    size_t comp2 = (comp + 1) % 2;

    if (!mBfield->useCompInMap(0, cell)) return false;
    cell[comp2]--;
    if (!mBfield->useCompInMap(0, cell)) return false;
    cell[comp2]++;
  }

  return true;
}

template<size_t DIM>
MxComplex MxYeeFitEField<DIM>::getCompFactor(size_t comp,
MxDimVector<int, DIM> cell) const {

  MxComplex res = MxYeeElecFieldBase<DIM>::getCompFactor(comp, cell);

  if (!this->useCompInMap(comp, cell)) return 0.0;
  else return res;
}


#if 0
template<size_t DIM>
void MxYeeFitEField<DIM>::setBCs(Teuchos::ParameterList bcList) {

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
MxDimVector<int, DIM> MxYeeFitEField<DIM>::getInteriorComp(size_t comp, MxDimVector<int, DIM> cell) const {
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
size_t MxYeeFitEField<DIM>::globCompIndx(size_t comp, const MxDimVector<int, DIM> & cell) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  return comp + this->numComps * this->grid->cellToGlobalIndx(newCell);
}

template<size_t DIM>
double MxYeeFitEField<DIM>::calcCompFrac(size_t comp, MxDimVector<int, DIM> cell, const MxShape<DIM> & aShape) const {
  MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
  MxDimVector<double, DIM> p(this->grid->nodeCoord(newCell) + this->compCellCoords[comp]);
  return this->compPtopes[comp]->volumeFraction(aShape, p);
}

template<size_t DIM>
MxComplex MxYeeFitEField<DIM>::getCompFactor(size_t comp, MxDimVector<int, DIM> cell) const {
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
bool MxYeeFitEField<DIM>::useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const {
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
template class MxYeeFitEField<1>;
template class MxYeeFitEField<2>;
template class MxYeeFitEField<3>;
