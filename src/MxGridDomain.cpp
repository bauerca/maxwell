
#include "MxGridDomain.h"
#include "MxUtil.hpp"
#include "MxDimVector.hpp"
#include "MxGrid.h"
#include "MxGridDomainIter.hpp"


template<size_t DIM>
MxGridDomain<DIM>::MxGridDomain(const MxGrid<DIM> * aGrid, const MxDimVector<int, DIM> & lowerBounds, const MxDimVector<int, DIM> & upperBounds, size_t guardCells) : 
grid(aGrid), lb(lowerBounds), ub(upperBounds), guard(guardCells), maxSize_t(-1) {
  interiorResolution = ub - lb;
  fullResolution = interiorResolution + MxDimVector<int, DIM>(2 * guard);
  numInteriorCells = interiorResolution.prod();
  numFullCells = fullResolution.prod();
  //numXCells = numInteriorCells - numMyCells;
  
}

template<size_t DIM>
MxDimVector<int, DIM> MxGridDomain<DIM>::interiorIndxToCell(int indx) const {
  MxDimVector<int, DIM> res;
  int factor = 1;
  for (size_t i = DIM - 1; i != maxSize_t; i--) {
    res[i] = (indx / factor) % interiorResolution[i] + lb[i];
    factor *= interiorResolution[i];
  }
  return res;
}


template<size_t DIM>
int MxGridDomain<DIM>::cellToInteriorIndx(const MxDimVector<int, DIM> & cell) const {
  int res = 0;
  int factor = 1;
  int gci, lbi, ubi;
  for (size_t i = DIM - 1; i != maxSize_t; i--) {
    gci = cell[i];
    lbi = lb[i];
    ubi = ub[i];

    // must make bounds check
    if ((gci < lbi) || (gci >= ubi)) {
      std::cout << "MxGridDomain: global cell input is out of bounds of interior cell block. Returning -1."; throw 1; }

    res += (gci - lbi) * factor;
    factor *= interiorResolution[i];
  }
  return res;
}



template<size_t DIM>
MxDimVector<int, DIM> MxGridDomain<DIM>::fullIndxToCell(int indx) const {
  MxDimVector<int, DIM> res;
  int factor = 1;
  for (size_t i = DIM - 1; i != maxSize_t; i--) {
    res[i] = (indx / factor) % fullResolution[i] + (lb[i] - guard);
    factor *= fullResolution[i];
  }
  return res;
}


template<size_t DIM>
int MxGridDomain<DIM>::cellToFullIndx(const MxDimVector<int, DIM> & cell) const {
  int res = 0;
  int factor = 1;
  int gci, lbi, ubi;
  for (size_t i = DIM - 1; i != maxSize_t; i--) {
    gci = cell[i];
    lbi = lb[i] - guard;
    ubi = ub[i] + guard;

    // must make bounds check
    if ((gci < lbi) || (gci >= ubi)) {
      std::cout << "MxGridDomain: global cell input is out of bounds of extended cell block\n";
      std::cout << "MxGridDomain: cell: "; cell.print();
      std::cout << "MxGridDomain: extended block lb: "; (lb - MxDimVector<int, DIM>(guard)).print();
      std::cout << "MxGridDomain: extended block ub: "; (ub + MxDimVector<int, DIM>(guard)).print();
      throw 1;
    }

    res += (gci - lbi) * factor;
    factor *= fullResolution[i];
  }
  return res;
}

template<size_t DIM>
MxGridDomainIter<DIM> MxGridDomain<DIM>::interiorBegin() const {
  MxGridDomainIter<DIM> res(grid);
  res.interiorDomIndx = 0;
  res.cell = this->interiorIndxToCell(res.interiorDomIndx);
  res.bump(0, 0); // initializes all other iter members
  return res;
}

template<size_t DIM>
MxGridDomainIter<DIM> MxGridDomain<DIM>::interiorEnd() const {
  MxGridDomainIter<DIM> res(grid);
  res.interiorDomIndx = numInteriorCells - 1;
  res.cell = this->interiorIndxToCell(res.interiorDomIndx);
  res.bump(0, 0); // initializes all other iter members
  return res;
}

template<size_t DIM>
MxGridDomainIter<DIM> MxGridDomain<DIM>::fullBegin() const {
  MxGridDomainIter<DIM> res(grid);
  res.fullDomIndx = 0;
  res.cell = this->fullIndxToCell(res.fullDomIndx);
  grid->interiorizeCell(res.cell);
  res.node = grid->nodeCoord(res.cell);
  res.globIndx = grid->cellToGlobalIndx(res.cell);
  return res;
}

template<size_t DIM>
MxGridDomainIter<DIM> MxGridDomain<DIM>::fullEnd() const {
  MxGridDomainIter<DIM> res(grid);
  res.fullDomIndx = numFullCells - 1;
  res.cell = this->fullIndxToCell(res.fullDomIndx);
  grid->interiorizeCell(res.cell);
  res.node = grid->nodeCoord(res.cell);
  res.globIndx = grid->cellToGlobalIndx(res.cell);
  return res;
}

// explicit instantiation

template class MxGridDomain<1>;
template class MxGridDomain<2>;
template class MxGridDomain<3>;
