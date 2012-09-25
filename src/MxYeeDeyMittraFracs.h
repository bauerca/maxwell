#ifndef MX_YEE_DM_FRACS
#define MX_YEE_DM_FRACS

#include "MxGridField.hpp"
#include "MxOperator.hpp"
#include "MxCrsMatrix.hpp"

#include "Teuchos_RCP.hpp"
#include "Epetra_CrsMatrix.h"

template<size_t DIM, typename Scalar>
class MxYeeDeyMittraFracs : public MxOperator<DIM>, public MxCrsMatrix<Scalar> {
  public:
    MxYeeDeyMittraFracs(const MxGridField<DIM> * aGridField,
      const MxGrid<DIM> * aGrid, bool getInverse = false,
      double aMinFrac = 0.0, bool randomize = false,
      double randomScale = 1.0);

    virtual ~MxYeeDeyMittraFracs() {};

    virtual const MxGridField<DIM> * getDomainField() const {return field;}

    virtual const MxGridField<DIM> * getRangeField() const {return field;}

    virtual const MxGrid<DIM> * getGrid() const {return grid;}

    //virtual const Epetra_CrsMatrix * getMatrixPtr() const {return matrix.get();}

    //virtual Epetra_CrsMatrix getMatrixCopy() const {return *matrix;}

    //virtual void clearMatrix() {matrix = Teuchos::null;}

  private:
    double minFrac;

    bool inverse;

    bool mRandomize;

    double mRandomScale;

    const MxGridField<DIM> * field;

    const MxGrid<DIM> * grid;

    void setMatrix();

    void setMatrixInverse();

    //Teuchos::RCP<Epetra_CrsMatrix> matrix;

};

#endif
