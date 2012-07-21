#ifndef MX_GRID_DOMAIN_ITER
#define MX_GRID_DOMAIN_ITER

#include "MxDimVector.hpp"
#include "MxGrid.h"
#include "MxGridDomain.h"


template<size_t DIM>
class MxGridDomainIter {
  public:
    MxGridDomainIter() {};

    MxGridDomainIter(const MxGrid<DIM> * aGrid) : grid(aGrid), domain(&aGrid->getDomain()) {};

    void fullBump();

    void interiorBump();

    void bump(size_t dir, int mag);

    bool operator==(const MxGridDomainIter<DIM> & iter) const {return (globIndx == iter.globIndx) and (fullDomIndx == iter.fullDomIndx);}
    bool operator!=(const MxGridDomainIter<DIM> & iter) const {return (globIndx != iter.globIndx) or (fullDomIndx != iter.fullDomIndx);}

    MxDimVector<int, DIM> getCell() const {return cell;}
    MxDimVector<double, DIM> getNode() const {return node;}
    size_t getGlobalIndx() const {return globIndx;}
    size_t getInteriorDomainIndx() const {return interiorDomIndx;}
    size_t getFullDomainIndx() const {return fullDomIndx;}

    void print() const;

    friend class MxGridDomain<DIM>;

  private:
    const MxGrid<DIM> * grid;
    const MxGridDomain<DIM> * domain;

    MxDimVector<double, DIM> node;
    MxDimVector<int, DIM> cell;
    size_t globIndx;
    size_t interiorDomIndx;
    size_t fullDomIndx;

};

template<size_t DIM>
inline
void MxGridDomainIter<DIM>::fullBump() {
  fullDomIndx++;
  cell = domain->fullIndxToCell(fullDomIndx);
  grid->interiorizeCell(cell);
  //interiorDomIndx = 0; // interiorDomIndx invalidated after fullBump
  globIndx = grid->cellToGlobalIndx(cell);
  node = grid->nodeCoord(cell);
}

template<size_t DIM>
inline
void MxGridDomainIter<DIM>::interiorBump() {
  interiorDomIndx++;
  cell = domain->interiorIndxToCell(interiorDomIndx);
  fullDomIndx = domain->cellToFullIndx(cell);
  globIndx = grid->cellToGlobalIndx(cell);
  node = grid->nodeCoord(cell);
}

template<size_t DIM>
inline
void MxGridDomainIter<DIM>::bump(size_t dir, int mag) {
  cell[dir] += mag;
  grid->interiorizeCell(cell);
  node = grid->nodeCoord(cell);
  fullDomIndx = domain->cellToFullIndx(cell);
  //interiorDomIndx = domain->cellToInteriorIndx(cell);
  globIndx = grid->cellToGlobalIndx(cell);
}

template<size_t DIM>
inline
void MxGridDomainIter<DIM>::print() const {
  std::cout << "MxGridDomainIter:\n"
            << "  Full Domain Index: " << fullDomIndx << "\n"
            << "  Interior Domain Index: " << interiorDomIndx << "\n"
            << "  Global Index: " << globIndx << "\n"
            << "  Grid Cell: "; cell.print();
  std::cout << "  Grid Node: "; node.print();
}


#endif
