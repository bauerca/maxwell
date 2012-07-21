#include "MxPML.h"

#include "MxGridFieldIter.hpp"
#include "MxEMSim.h"
#include "MxIO.h"

#include "Epetra_MultiVector.h"


template<size_t DIM>
void MxPML<DIM>::alphaMap(MxEMSim<DIM> const * sim, MxGridField<DIM> const * field) {
  MxMultiVector<double> mv(field->getMap(), 1);

  MxGridFieldIter<DIM> iter(field);

  int row;
  MxDimVector<int, DIM> cell;
  MxDimVector<double, DIM> coord;

  MxPML<DIM> const * pml;

  double alpha;

  for (iter.begin(); !iter.atEnd(); iter.bump()) {
    row = iter.getGlobCompIndx();
    cell = iter.getCell();
    coord = iter.getCoord();

    alpha = 0;
    for (size_t i = 0; i < sim->numPMLs(); ++i) {
      pml = sim->getPML(i);
      if (pml->getShape()->func(coord) >= 0)
        alpha += pml->getAlpha(coord);
    }

    mv.replaceGlobalValue(row, 0, alpha);
  }

  MxIO<DIM, double> io(sim->getGrid().getComm());
  io.save(mv, *field, "alphaMap");
}


template class MxPML<1>;
template class MxPML<2>;
template class MxPML<3>;
