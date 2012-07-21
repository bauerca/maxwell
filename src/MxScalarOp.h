#ifndef MX_SCALAR_OP
#define MX_SCALAR_OP

#include "MxTypes.h"
#include "MxOperator.hpp"
#include "MxGridField.hpp"

#include "Epetra_CrsMatrix.h"

template<size_t DIM>
class MxScalarOp : public MxOperator<DIM>, public Epetra_CrsMatrix {
  public:
    MxScalarOp(MxComplex value, MxGridField<DIM> const * field, bool imag = false);

    virtual ~MxScalarOp() {};

  private:

    MxComplex mValue;

    MxGridField<DIM> const * mField;

    void build();

};

#endif
