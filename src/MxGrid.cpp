#include <iostream>
#include <cmath>

#include "MxGrid.h"
#include "MxGridDomainIter.hpp"
#include "MxUtil.hpp"
#include "MxMap.hpp"


template<size_t DIM>
MxGrid<DIM>::MxGrid(const MxDimVector<double, DIM> & gridOrigin, 
const MxDimVector<int, DIM> & gridResolution, 
const MxDimVector<double, DIM> & gridSize, 
RCP<MxComm> comm) : 
resolution(gridResolution), size(gridSize), origin(gridOrigin), numCells(gridResolution.prod()), 
mComm(comm), maxSize_t(-1) {
  // set the cell size and bounds
  for (size_t i = 0; i < DIM; i++) {
    cellSize[i] = size[i] / double(resolution[i]);
    ub[i] = origin[i] + size[i];
    lb[i] = origin[i];
  }
  this->partition();
}


template<size_t DIM>
MxGrid<DIM>::MxGrid(Teuchos::XMLObject const & node,
RCP<MxComm> comm) :
resolution(0), size(0), origin(0), numCells(0), 
mComm(comm), maxSize_t(-1) {
  MxDimVector<int, DIM> n;
  MxDimVector<double, DIM> l, o;

  resolution.strFill(MxUtil::XML::getAttr("resolution", node, "[16, 16, 16]"));
  size.strFill(MxUtil::XML::getAttr("dimensions", node, "[1.0, 1.0, 1.0]"));
  origin.strFill(MxUtil::XML::getAttr("origin", node, "[0.0, 0.0, 0.0]"));

  numCells = resolution.prod();

  for (size_t i = 0; i < DIM; ++i) {
    cellSize[i] = size[i] / double(resolution[i]);
    ub[i] = origin[i] + size[i];
    lb[i] = origin[i];
  }
  this->partition();
}


template<size_t DIM>
MxMap MxGrid<DIM>::fieldMap(int numComps) const {
  MxGridDomain<DIM> domain = getGridDomain(0);

#if 1
  std::vector<MxIndex> globCompInds;
  globCompInds.reserve(numComps * domain.getNumInteriorCells());

  // using a grid domain iterator
  MxGridDomainIter<DIM> gdIter = domain.interiorBegin();
  MxGridDomainIter<DIM> gdIterEnd = domain.interiorEnd();
  MxDimVector<int, DIM> cell;
  for (gdIter; gdIter != gdIterEnd; gdIter.interiorBump()) {
    cell = gdIter.getCell();
    for (int comp = 0; comp < numComps; ++comp) {
      globCompInds.push_back(comp + numComps * this->cellToGlobalIndx(cell));
    }
  }
#endif

#if 0
  int numIntCells = domain->getNumInteriorCells();
  int numIntComps = numComps * numIntCells;

  std::vector<int> globCompInds;
  globCompInds.reserve(numIntComps);

  MxDimVector<int, DIM> cell;
  int globCellIndx, globCompIndx;
  for (int i = 0; i < numIntCells; i++) {
    cell = domain->interiorIndxToCell(i);
    globCellIndx = grid->cellToGlobalIndx(cell);
    for (int comp = 0; comp < numComps; ++comp) {
      if (this->useCompInMap(comp, cell)) {
        globCompIndx = comp + numComps * globCellIndx;
        globCompInds.push_back(globCompIndx);
      }
    }
    //addGlobalCompIndices(cell, gridRes, *domain, globCompInds);
  }
#endif

  size_t numGlobElems = numComps * (MxDimVector<int, DIM>(1) + resolution).prod();
  //return MxMap(numGlobElems, globCompInds.size(), &globCompInds[0], 0, *mComm);
  return MxMap(numGlobElems, globCompInds, mComm);
}



template<size_t DIM>
void MxGrid<DIM>::partition() {
  // Domain decomposition
  size_t Np = mComm->numProc();
  //Np = 16 * 4;
  if (mComm->myPID() == 0) std::cout << "Number of processes: " << Np << "\n";

  std::vector<std::vector<size_t> > allFacs, facs, perms, decomps;
  //facs = MxUtil::factorizations(Np, DIM);
  allFacs = MxUtil::unorderedFactorizations(Np);
  size_t numFacs = allFacs.size();
  //for (int i = 0; i < numFacs; ++i)
  //  MxUtil::printStdVector(allFacs[i]);

  // add 1's to factorizations with fewer than DIM terms and drop
  // those with more than DIM.
  for (int i = 0; i < numFacs; ++i) {
    size_t terms = allFacs[i].size();
    if (terms < DIM) {
      std::vector<size_t> fac(allFacs[i]);
      for (int j = 0; j < DIM - terms; ++j)
        fac.push_back(1);
      facs.push_back(fac);
    }
    else if (terms == DIM)
      facs.push_back(allFacs[i]);
  }
  numFacs = facs.size();
  //for (int i = 0; i < numFacs; ++i)
  //  MxUtil::printStdVector(facs[i]);
  
  // get all permutations of each factorization of the number of processes
  if (mComm->myPID() == 0) std::cout << "perms of facs:\n";
  for (size_t i = 0; i < numFacs; i++) {
    perms = MxUtil::permutations<size_t>(facs[i]);
    for (size_t j = 0; j < perms.size(); j++) {
      decomps.push_back(perms[j]);
      if (mComm->myPID() == 0) MxUtil::printStdVector(perms[j]);
    }
  }
  if (mComm->myPID() == 0) std::cout << "end perms of facs.\n\n";

  size_t numDecomps = decomps.size();

  // for each decomposition, test the surface-to-volume
  // ratio of the resulting domain size.
  double ratio = 0, minRatio = 0; 
  size_t idealFacIndx = 0;
  for (size_t i = 0; i < numDecomps; i++) {
    // get all permutations of current decomposition
    //if (mComm->myPID() == 0) MxUtil::printStdVector(decomps[i]);
    ratio = 0;
    for (size_t j = 0; j < DIM; j++)
      ratio += double(decomps[i][j]) / double(resolution[j]);
    //if (mComm->myPID() == 0) std::cout << "  surf/vol ratio: " << ratio << "\n";

    if ((i == 0) || (ratio < minRatio)) {
      minRatio = ratio;
      idealFacIndx = i;
    }
  }

  MxDimVector<int, DIM> procGrid; // grid of processes
  for (size_t i = 0; i < DIM; i++)
    procGrid[i] = decomps[idealFacIndx][i];

  domainDecomp = procGrid;
  // now create grid subdomain for this processor
  this->setGridDomain(procGrid);
}

template<size_t DIM>
MxGridDomain<DIM> MxGrid<DIM>::getGridDomain(int numGuard) const {
  int pid = mComm->myPID();

  MxDimVector<int, DIM> lb, ub;

  int factor = 1;
  int stdRes, xRes, proci;
  for (size_t i = DIM - 1; i != maxSize_t; i--) {
    stdRes = int(floor(double(resolution[i]) / double(domainDecomp[i]) + 0.5));
    xRes = resolution[i] - (domainDecomp[i] - 1) * stdRes;

    proci = (pid / factor) % domainDecomp[i];
    lb[i] = proci * stdRes;
    if (proci == domainDecomp[i] - 1)
      ub[i] = resolution[i] + 1;
    else
      ub[i] = lb[i] + stdRes;

    factor *= domainDecomp[i];
  }

  return MxGridDomain<DIM>(this, lb, ub, numGuard);
}

template<size_t DIM>
void MxGrid<DIM>::setGridDomain(const MxDimVector<int, DIM> & dd) {
  int pid = mComm->myPID();

  MxDimVector<int, DIM> lb, ub;

  int factor = 1;
  int stdRes, xRes, proci;
  for (size_t i = DIM - 1; i != maxSize_t; i--) {
    stdRes = int(floor(double(resolution[i]) / double(dd[i]) + 0.5));
    xRes = resolution[i] - (dd[i] - 1) * stdRes;

    proci = (pid / factor) % dd[i];
    lb[i] = proci * stdRes;
    if (proci == dd[i] - 1)
      ub[i] = resolution[i] + 1;
    else
      ub[i] = lb[i] + stdRes;

    factor *= dd[i];
  }

  myDomain = Teuchos::rcp(new MxGridDomain<DIM>(this, lb, ub));
}

template<size_t DIM>
void MxGrid<DIM>::interiorizeCell(MxDimVector<int, DIM> & cell) const {
  int gci;
  int ni;
  for (size_t i = 0; i < DIM; i++) {
    gci = cell[i];
    ni = resolution[i];
    if (gci >= ni)
      cell[i] = gci - ni;
    else if (gci < 0)
      cell[i] = ni + gci;
  }
}

template<size_t DIM>
void MxGrid<DIM>::print() const {
  if (mComm->myPID() == 0) {
    std::cout << "============= Grid =============\n"
              << "  Resolution: ";
    for (size_t i = 0; i < DIM; i++) {
      if (i == DIM - 1) std::cout << resolution[i] << "\n";
      else std::cout << resolution[i] << " x ";
    }
    std::cout << "  Origin: (";
    for (size_t i = 0; i < DIM; i++) {
      if (i == DIM - 1) std::cout << origin[i] << ")\n";
      else std::cout << origin[i] << ", ";
    }
    std::cout << "  Size: ";
    for (size_t i = 0; i < DIM; i++) {
      if (i == DIM - 1) std::cout << size[i] << "\n";
      else std::cout << size[i] << " x ";
    }
    std::cout << "  Cell Size: ";
    for (size_t i = 0; i < DIM; i++) {
      if (i == DIM - 1) std::cout << cellSize[i] << "\n";
      else std::cout << cellSize[i] << " x ";
    }
    std::cout << "  Domain decomposition: ";
    for (size_t i = 0; i < DIM; i++) {
      if (i == DIM - 1) std::cout << domainDecomp[i] << "\n";
      else std::cout << domainDecomp[i] << " x ";
    }
    std::cout << "================================\n";
  }
}

template<size_t DIM>
MxDimVector<double, DIM> MxGrid<DIM>::nodeCoord(const MxDimVector<int, DIM> & cell) const {
  MxDimVector<double, DIM> res;
  for (size_t i = 0; i < DIM; i++)
    res[i] = origin[i] + double(cell[i]) * cellSize[i];
  return res;
}

template<>
std::vector<MxDimVector<int, 1> > MxGrid<1>::getCellBlock(MxDimVector<int, 1> originCell, int size) const {
  std::vector<MxDimVector<int, 1> > cellBlock(size, originCell);
  int cellIndx = 0;
  for (int i = 0; i < size; i++) {
    cellBlock[cellIndx][0] += i;
    cellIndx++;
  }
  return cellBlock;
}

template<>
std::vector<MxDimVector<int, 2> > MxGrid<2>::getCellBlock(MxDimVector<int, 2> originCell, int size) const {
  std::vector<MxDimVector<int, 2> > cellBlock(size * size, originCell);
  int cellIndx = 0;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      cellBlock[cellIndx][0] += i;
      cellBlock[cellIndx][1] += j;
      cellIndx++;
    }
  }
  return cellBlock;
}

template<>
std::vector<MxDimVector<int, 3> > MxGrid<3>::getCellBlock(MxDimVector<int, 3> originCell, int size) const {
  std::vector<MxDimVector<int, 3> > cellBlock(size * size * size, originCell);
  int cellIndx = 0;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      for (int k = 0; k < size; k++) {
        cellBlock[cellIndx][0] += i;
        cellBlock[cellIndx][1] += j;
        cellBlock[cellIndx][2] += k;
        cellIndx++;
      }
    }
  }
  return cellBlock;
}

template<size_t DIM>
bool MxGrid<DIM>::operator==(MxGrid<DIM> const & theGrid) const {
  for (size_t i = 0; i < DIM; ++i) {
    if (resolution[i] != theGrid.resolution[i])
      return false;
    if (size[i] != theGrid.size[i])
      return false;
    if (origin[i] != theGrid.origin[i])
      return false;
  }
  return true;
}

template class MxGrid<1>;
template class MxGrid<2>;
template class MxGrid<3>;

