#include "MxScalarOp.h"
#include "MxGridFieldIter.hpp"


template<size_t DIM>
MxScalarOp<DIM>::MxScalarOp(MxComplex value, MxGridField<DIM> const * field, bool imag) :
Epetra_CrsMatrix(Copy, field->getMap(), 1, true), mValue(value), mField(field) {
  this->imaginary = imag;
  build();
}


template<size_t DIM>
void MxScalarOp<DIM>::build() {
  int row;
  double val;

  MxGridFieldIter<DIM> iter(field);
  for (iter.begin(); !iter.atEnd(); iter.bump()) {
    row = iter.getGlobCompIndx();
    if (this->imaginary)
      val = mValue.imag();
    else
      val = mValue.real();
    Epetra_CrsMatrix::InsertGlobalValues(row, 1, &val, &row);
  }
  Epetra_CrsMatrix::FillComplete(mField->getMap(), mField->getMap());
}


template class MxScalarOp<1>;
template class MxScalarOp<2>;
template class MxScalarOp<3>;
