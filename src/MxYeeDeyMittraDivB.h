#ifndef MX_YEE_DM_DIV_B
#define MX_YEE_DM_DIV_B

#include "MxGrid.h"
#include "MxGridField.hpp"
#include "MxEMSim.h"
#include "MxOperator.hpp"
#include "MxCrsMatrix.hpp"

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

template<size_t DIM, typename Scalar>
class MxYeeDeyMittraDivB : public MxOperator<DIM>, public MxCrsMatrix<Scalar> {
  public:
    MxYeeDeyMittraDivB(RCP<MxEMSim<DIM> > theSim);

    virtual const MxGridField<DIM> * getDomainField() const {return bfield;}

    virtual const MxGridField<DIM> * getRangeField() const {return psifield;}

    virtual const MxGrid<DIM> * getGrid() const {return grid;}

  private:
    const MxGridField<DIM> * bfield;

    const MxGridField<DIM> * psifield;

    const MxGrid<DIM> * grid;

    void setMatrix();

};

#endif
