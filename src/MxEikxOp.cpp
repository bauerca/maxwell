#include "MxEikxOp.h"


template<size_t DIM>
MxEikxOp<DIM>::MxEikxOp(MxDimVector<double, DIM> blochK,
double sign, double shift,
const MxGridField<DIM> * aGridField,
const MxGrid<DIM> * aGrid) :
MxGridFieldFuncOp<DIM, MxComplex>(aGridField, aGrid),
I(MxComplex(0.0, 1.0)), mBlochK(blochK), mSign(sign), mShift(shift) {
  this->build();
}


template<size_t DIM>
MxComplex MxEikxOp<DIM>::func(size_t comp,
MxDimVector<int, DIM> cell) const {
  MxDimVector<double, DIM> coord(this->field->getCompCoord(comp, cell));
  return exp(mSign * I * (coord.dot(mBlochK) + mShift));
}

template class MxEikxOp<1>;
template class MxEikxOp<2>;
template class MxEikxOp<3>;
