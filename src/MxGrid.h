#ifndef MX_GRID
#define MX_GRID

#include "MxComm.hpp"
#include "MxDimVector.hpp"
#include "MxGridDomain.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLObject.hpp"

class MxMap;

template<size_t DIM>
class MxGrid {
  public:
    MxGrid(const MxDimVector<double, DIM> & gridOrigin,
      const MxDimVector<int, DIM> & gridResolution,
      const MxDimVector<double, DIM> & gridSize,
      RCP<MxComm> comm);

    MxGrid(Teuchos::XMLObject const & node,
      RCP<MxComm> comm);

    //int cellToDomainInteriorIndx(const MxDimVector<int, DIM> & cell) const;

    //int cellToDomainFullIndx(const MxDimVector<int, DIM> & cell) const;

    int cellToGlobalIndx(const MxDimVector<int, DIM> & cell) const;

    MxDimVector<int, DIM> globalIndxToCell(int indx) const;

    void interiorizeCell(MxDimVector<int, DIM> & cell) const;

    void print() const;

    size_t getNumCells() const {return numCells;}

    MxDimVector<double, DIM> getCellSize() const {return cellSize;}

    MxDimVector<int, DIM> getResolution() const {return resolution;}

    MxDimVector<double, DIM> getSize() const {return size;}

    MxDimVector<double, DIM> getOrigin() const {return origin;}

    MxDimVector<double, DIM> getUpperBounds() const {return ub;}

    MxDimVector<double, DIM> getLowerBounds() const {return lb;}

    const MxGridDomain<DIM> & getDomain() const {return *myDomain;}

    MxGridDomain<DIM> getGridDomain(int numGuard = 1) const;

    std::vector<MxDimVector<int, DIM> > getCellBlock(MxDimVector<int, DIM> originCell, int size) const;

    MxDimVector<double, DIM> nodeCoord(const MxDimVector<int, DIM> & cell) const;
    
    MxMap fieldMap(int numComps) const;

    //int globalIndxToDomainInteriorIndx(int indx) const;

    RCP<MxComm> getComm() const {return mComm;}

    bool operator==(MxGrid<DIM> const & theGrid) const;
    
    friend class MxGridDomainIter<DIM>;


  private:
    MxDimVector<int, DIM> resolution;

    MxDimVector<double, DIM> size;

    MxDimVector<double, DIM> origin;

    MxDimVector<double, DIM> cellSize;

    MxDimVector<double, DIM> ub, lb;

    MxDimVector<int, DIM> domainDecomp;

    Teuchos::RCP<MxGridDomain<DIM> > myDomain;

    size_t numCells;

    void partition();

    void setGridDomain(const MxDimVector<int, DIM> & dd);

    RCP<MxComm> mComm;

    size_t maxSize_t;

};

template<size_t DIM>
inline
int MxGrid<DIM>::cellToGlobalIndx(const MxDimVector<int, DIM> & cell) const {
  int res = 0;
  int factor = 1;
  int gci, ni; 
  for (size_t i = DIM - 1; i != maxSize_t; i--) {
    gci = cell[i];
    ni = resolution[i] + 1;
    if (gci >= ni)
      res += (gci - ni) * factor; // always periodic
    else if (gci < 0)
      res += (ni + gci) * factor; // always periodic
    else
      res += gci * factor;
    factor *= ni;
  }
  return res;
}

template<size_t DIM>
inline
MxDimVector<int, DIM> MxGrid<DIM>::globalIndxToCell(int indx) const {
  MxDimVector<int, DIM> res;
  int factor = 1;
  int ni;
  for (size_t i = DIM - 1; i != maxSize_t; i--) {
    ni = resolution[i] + 1;
    res[i] = (indx / factor) % ni;
    factor *= ni;
  }
  return res;
}


#endif
