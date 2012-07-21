#ifndef MX_GRID_FIELD_FUNC_OP
#define MX_GRID_FIELD_FUNC_OP

#include "MxTypes.h"

#include "MxGridField.hpp"
#include "MxOperator.hpp"
#include "MxCrsMatrix.hpp"

#include "Teuchos_RCP.hpp"

template<size_t DIM, typename Scalar>
class MxGridFieldFuncOp : public MxOperator<DIM>, public MxCrsMatrix<Scalar> {
  public:
    MxGridFieldFuncOp(const MxGridField<DIM> * aGridField,
        const MxGrid<DIM> * aGrid);

    virtual ~MxGridFieldFuncOp() {};

    virtual const MxGridField<DIM> * getDomainField() const {return field;}

    virtual const MxGridField<DIM> * getRangeField() const {return field;}

    virtual const MxGrid<DIM> * getGrid() const {return grid;}

  protected:

    virtual Scalar func(size_t comp, MxDimVector<int, DIM> cell) const {
      return ScalarTraits<Scalar>::one();
    }

    virtual void build();

    //virtual const Epetra_CrsMatrix * getMatrixPtr() const {return matrix.get();}

    //virtual Epetra_CrsMatrix getMatrixCopy() const {return *matrix;}

    //virtual void clearMatrix() {matrix = Teuchos::null;}

    const MxGridField<DIM> * field;

    const MxGrid<DIM> * grid;


};

#endif
