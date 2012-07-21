#ifndef MX_EM_SIM_HIERARCHY
#define MX_EM_SIM_HIERARCHY

#include <vector>

#include "MxGrid.h"
#include "MxEMSim.h"

#include "Teuchos_RCP.hpp"

template<size_t DIM>
class MxEMSimHierarchy {
  public:
    explicit MxEMSimHierarchy(const MxEMSim<DIM> * aSim, int approxLevels = 1);

    const MxEMSim<DIM> * getSim(int level) const;

    int getNumLevels() const {return levels;}


  private:
    int levels;

    const MxEMSim<DIM> * fineSim;

    std::vector<Teuchos::RCP<MxGrid<DIM> > > coarseGrids;

    std::vector<Teuchos::RCP<MxEMSim<DIM> > > coarseSims;

    void coarsen(const MxEMSim<DIM> & aSim, double factor, MxEMSim<DIM> *& simPtr, MxGrid<DIM> *& gridPtr);


};

#endif
