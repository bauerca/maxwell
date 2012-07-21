#ifndef MX_GRID_DOMAIN
#define MX_GRID_DOMAIN

//template<typename T, size_t DIM> MxDimVector;
#include "MxDimVector.hpp"

template<size_t DIM> class MxGridDomainIter;
template<size_t DIM> class MxGrid;

template<size_t DIM>
class MxGridDomain {
  public:
    MxGridDomain(const MxGrid<DIM> * aGrid, const MxDimVector<int, DIM> & lowerBounds, const MxDimVector<int, DIM> & upperBounds, size_t guardCells = 1);

    MxDimVector<int, DIM> interiorIndxToCell(int indx) const;

    int cellToInteriorIndx(const MxDimVector<int, DIM> & cell) const;

    MxDimVector<int, DIM> fullIndxToCell(int indx) const;

    int cellToFullIndx(const MxDimVector<int, DIM> & cell) const;

    size_t getNumFullCells() const {return numFullCells;}

    size_t getNumInteriorCells() const {return numInteriorCells;}

    MxGridDomainIter<DIM> interiorBegin() const;

    MxGridDomainIter<DIM> interiorEnd() const;

    MxGridDomainIter<DIM> fullBegin() const;

    MxGridDomainIter<DIM> fullEnd() const;

    MxDimVector<int, DIM> getInteriorResolution() const {return interiorResolution;}

    MxDimVector<int, DIM> getFullResolution() const {return fullResolution;}

    MxDimVector<int, DIM> getLowerBoundCell() const {return lb;}

    MxDimVector<int, DIM> getUpperBoundCell() const {return ub;}

    friend class MxGridDomainIter<DIM>;
    friend class MxGrid<DIM>;

  private:
    const MxGrid<DIM> * grid;

    // cell indices marking domain box boundaries
    MxDimVector<int, DIM> lb, ub;

    size_t guard;

    // number of cells in each dimension for interior block
    MxDimVector<int, DIM> interiorResolution;

    // number of cells in each dimension for entire block
    // including guard cells
    MxDimVector<int, DIM> fullResolution;

    size_t numInteriorCells;

    size_t numFullCells;

    size_t maxSize_t;
};



#endif
