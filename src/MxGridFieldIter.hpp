#ifndef MX_GRID_FIELD_ITER
#define MX_GRID_FIELD_ITER

#include "MxDimVector.hpp"
#include "MxGridField.hpp"


template<size_t DIM>
class MxGridFieldIter {
  public:
    MxGridFieldIter() {};

    MxGridFieldIter(const MxGridField<DIM> * aField) : 
    field(aField), localIndx(0), 
    numInds(aField->getMap()->getNodeNumIndices()), 
    globMapInds(aField->getMap()->getNodeIndexList()),
    cellSize(aField->grid->getCellSize()) {};

    inline void begin() {localIndx = 0;}

    inline void bump() {localIndx++;}

    inline bool atEnd() {return localIndx == numInds;}

    bool operator==(const MxGridFieldIter<DIM> & iter) const {return localIndx == iter.localIndx;}
    bool operator!=(const MxGridFieldIter<DIM> & iter) const {return localIndx != iter.localIndx;}

    MxDimVector<int, DIM> getCell();

    inline size_t getComp() {return globMapInds[localIndx] % field->getNumComps();}

    inline int getGlobCompIndx() {return globMapInds[localIndx];}

    inline int getLocalCompIndx() {return localIndx;}

    MxDimVector<double, DIM> getCoord();

    void print() const;

  private:
    const MxGridField<DIM> * field;

    int localIndx;

    int numInds;
    MxIndex const * globMapInds;

    MxDimVector<double, DIM> cellSize;
    MxDimVector<int, DIM> cell;
    MxDimVector<double, DIM> coord;
    size_t globCellIndx;

};

template<size_t DIM>
inline
MxDimVector<int, DIM> MxGridFieldIter<DIM>::getCell() {
  globCellIndx = globMapInds[localIndx] / field->getNumComps();
  cell = field->grid->globalIndxToCell(globCellIndx);
  return cell;
}


template<size_t DIM>
inline
void MxGridFieldIter<DIM>::print() const {
  std::cout << "MxGridFieldIter:\n"
            << "  Global Cell Index: " << globCellIndx << "\n"
            << "  Grid Cell: "; cell.print();
}

template<size_t DIM>
inline
MxDimVector<double, DIM> MxGridFieldIter<DIM>::getCoord() {
  cell = this->getCell();
  coord = field->grid->nodeCoord(cell) + field->compCellCoords[this->getComp()];
  return coord;
}

#endif
