
#include "MxEMSimHierarchy.h"

#include "MxDimVector.hpp"
#include "MxGrid.h"

template<size_t DIM>
MxEMSimHierarchy<DIM>::MxEMSimHierarchy(const MxEMSim<DIM> * aSim, int approxLevels) :
levels(approxLevels), fineSim(aSim) {

  MxEMSim<DIM> * simPtr;
  MxGrid<DIM> * gridPtr;
  if (approxLevels > 1) {
    // coarsen from fine simulation first
    coarsen(*fineSim, 2, simPtr, gridPtr);

    if (simPtr != 0) {
      coarseSims.push_back(Teuchos::rcp(simPtr));
      coarseGrids.push_back(Teuchos::rcp(gridPtr));
    }

    // coarsen from coarse simulations
    for (int i = 1; i < approxLevels - 1; ++i) {
      coarsen(*coarseSims[i - 1], 2, simPtr, gridPtr);

      if (simPtr != 0) {
        coarseSims.push_back(Teuchos::rcp(simPtr));
        coarseGrids.push_back(Teuchos::rcp(gridPtr));
      }
    }
    
  }
}


template<size_t DIM>
void MxEMSimHierarchy<DIM>::coarsen(const MxEMSim<DIM> & aSim, double factor, MxEMSim<DIM> *& simPtr, MxGrid<DIM> *& gridPtr) {

  const MxGrid<DIM> & grid = aSim.getGrid();

  MxDimVector<int, DIM> resolution;
  MxDimVector<double, DIM> size, origin, lb, ub;

  // halve the resolution (or approximately halve the resolution)
  resolution = grid.getResolution() / 2;

  const MxShape<DIM> * pec = aSim.getPEC();
  if (pec == 0) {
    origin = grid.getOrigin();
    size = grid.getSize();
  }
  else {
    lb = pec->boundingBox().first;
    ub = pec->boundingBox().second;

    // make coarse grid larger than conducting bounding box by 4 fine cells
    // (coarse grid will be larger than bbox by 2 cells if fine resolution is multiple
    // of two in each dimension)
    origin = lb - 4.0 * grid.getCellSize();
    size = (ub - lb) + 2.0 * 4.0 * grid.getCellSize();
  }

  //simPtr = new MxEMSim<DIM>(aSim);

  gridPtr = new MxGrid<DIM>(origin, resolution, size, &grid.getComm());
  gridPtr->print();

  simPtr = new MxEMSim<DIM>(gridPtr, aSim.getParameters());
  simPtr->setPEC(aSim.getPEC());
  for (size_t i = 0; i < aSim.numDielectrics(); ++i)
    simPtr->addDielectric(aSim.getDielectric(i));

  //simPtr->setGrid(gridPtr);
  simPtr->setup();
  //simPtr = new MxEMSim<DIM>(*gridPtr, aSim.getParameters());

}

template<size_t DIM>
const MxEMSim<DIM> * MxEMSimHierarchy<DIM>::getSim(int level) const {
  if (level == 0)
    return fineSim;
  else
    return coarseSims[level - 1].get();
}



template class MxEMSimHierarchy<1>;
template class MxEMSimHierarchy<2>;
template class MxEMSimHierarchy<3>;
