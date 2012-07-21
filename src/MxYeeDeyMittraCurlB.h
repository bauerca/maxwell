#ifndef MX_YEE_DM_CURL_B
#define MX_YEE_DM_CURL_B

#include "MxGrid.h"
#include "MxGridField.hpp"
#include "MxEMSim.h"
#include "MxOperator.hpp"
#include "MxCrsMatrix.hpp"

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

template<size_t DIM, typename Scalar>
class MxYeeDeyMittraCurlB : public MxOperator<DIM>, public MxCrsMatrix<Scalar> {
  public:
    MxYeeDeyMittraCurlB(RCP<MxEMSim<DIM> > theSim);
    
    virtual const MxGridField<DIM> * getDomainField() const {return efield;}

    virtual const MxGridField<DIM> * getRangeField() const {return bfield;}

    virtual const MxGrid<DIM> * getGrid() const {return grid;}


    //virtual void setMatrix();

    //virtual const Epetra_CrsMatrix * getMatrixPtr() const {return matrix.get();}

    //virtual Epetra_CrsMatrix getMatrixCopy() const {return *matrix;}

    //virtual void clearMatrix() {matrix = Teuchos::null;}

  private:
    const MxGridField<DIM> * efield;

    const MxGridField<DIM> * bfield;

    const MxGrid<DIM> * grid;

    void setMatrix2d();

    void setMatrix3d();

    //Teuchos::RCP<Epetra_CrsMatrix> matrix;
};

#endif
