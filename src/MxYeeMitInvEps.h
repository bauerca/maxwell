#ifndef MX_YEE_MIT_INV_EPS
#define MX_YEE_MIT_INV_EPS

#include <vector>

#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxDynMatrix.hpp"
#include "MxGridField.hpp"
#include "MxEMSim.h"
#include "MxDielectric.hpp"
#include "MxOperator.hpp"

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

template<size_t DIM>
class MxYeeMitInvEps : public MxOperator<DIM>, public Epetra_CrsMatrix {
  public:
    explicit MxYeeMitInvEps(const MxEMSim<DIM> * theSim, bool imag = false);
    
    virtual const MxGridField<DIM> * getDomainField() const {return efield;}

    virtual const MxGridField<DIM> * getRangeField() const {return efield;}

    virtual const MxGrid<DIM> * getGrid() const {return grid;}

    //virtual void setMatrix();

    //virtual const Epetra_CrsMatrix * getMatrixPtr() const {return matrix.get();}

    //virtual Epetra_CrsMatrix getMatrixCopy() const {return *matrix;}

    //virtual void clearMatrix() {matrix = Teuchos::null;}

    Epetra_CrsMatrix cellAveInvEpsOperator() const;

  private:

    const MxEMSim<DIM> * sim;

    const MxGridField<DIM> * efield;

    const MxGridField<DIM> * dfield;

    const MxDielectric<DIM> * diel;

    const MxGrid<DIM> * grid;

    bool epsIsDiagonal(const MxDimMatrix<double, 3> & epsTensor) const;

    std::vector<MxDimMatrix<double, DIM> > epsOff;

    std::vector<std::vector<MxDimMatrix<double, DIM> > > epsDiag;

    void setMatrix2dTM();

    void setMatrix2dTE3d();

    void fillAveEps();

    MxPolytope<DIM> * getVoxel() const;

    MxDimMatrix<double, DIM> aveInvEps(MxDimVector<double, DIM> pt, double vfrac, MxDimMatrix<double, DIM> epsIn, MxDimMatrix<double, DIM> epsOut) const;

    void update(int comp0, MxDimVector<int, DIM> cell, std::vector<double> & vals, std::vector<int> & inds, std::vector<double> & factors);

    //Teuchos::RCP<Epetra_CrsMatrix> matrix;
};

#endif
