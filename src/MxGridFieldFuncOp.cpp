
#include "MxGridFieldFuncOp.h"
#include "MxDimVector.hpp"
#include "MxGridFieldIter.hpp"



template<size_t DIM, typename Scalar>
MxGridFieldFuncOp<DIM, Scalar>::MxGridFieldFuncOp(
const MxGridField<DIM> * aGridField, 
const MxGrid<DIM> * aGrid) : 
MxCrsMatrix<Scalar>(aGridField->getMap(), aGridField->getMap()),
field(aGridField), grid(aGrid) {}


template<size_t DIM, typename Scalar>
void MxGridFieldFuncOp<DIM, Scalar>::build() {
  MxGridFieldIter<DIM> fieldIter(field);

  double minNonzeroFrac = 1.0;
  int numZeroFracs = 0;

  MxIndex ind;
  size_t comp;
  Scalar val;
  MxDimVector<int, DIM> cell;
  for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
    ind = fieldIter.getGlobCompIndx();
    comp = fieldIter.getComp();
    cell = fieldIter.getCell();
    val = this->func(comp, cell);
    MxCrsMatrix<Scalar>::insertRowValues(ind, 1, &ind, &val);
  }
  MxCrsMatrix<Scalar>::fillComplete(field->getMap(), field->getMap());
}


template class MxGridFieldFuncOp<1, double>;
template class MxGridFieldFuncOp<2, double>;
template class MxGridFieldFuncOp<3, double>;
template class MxGridFieldFuncOp<1, MxComplex>;
template class MxGridFieldFuncOp<2, MxComplex>;
template class MxGridFieldFuncOp<3, MxComplex>;
