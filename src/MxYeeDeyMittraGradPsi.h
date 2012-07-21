#ifndef MX_YEE_DM_GRAD_PSI
#define MX_YEE_DM_GRAD_PSI

#include "MxGrid.h"
#include "MxGridField.hpp"
#include "MxEMSim.h"
#include "MxOperator.hpp"
#include "MxCrsMatrix.hpp"

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

template<size_t DIM, typename Scalar>
class MxYeeDeyMittraGradPsi : public MxOperator<DIM>, public MxCrsMatrix<Scalar> {
  public:
    MxYeeDeyMittraGradPsi(RCP<MxEMSim<DIM> > theSim);

    virtual const MxGridField<DIM> * getDomainField() const {return psifield;}

    virtual const MxGridField<DIM> * getRangeField() const {return bfield;}

    virtual const MxGrid<DIM> * getGrid() const {return grid;}

  private:
    const MxGridField<DIM> * psifield;

    const MxGridField<DIM> * bfield;

    const MxGrid<DIM> * grid;

    void setMatrix();

};

#endif
