#ifndef MX_YEE_FIT_MU
#define MX_YEE_FIT_MU

#include <vector>

#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxDynMatrix.hpp"
#include "MxGridField.hpp"
#include "MxEMSim.h"
#include "MxDielectric.hpp"
#include "MxOperator.hpp"
#include "MxCrsMatrix.hpp"

template<size_t DIM, typename Scalar>
class MxYeeFitMu : public MxOperator<DIM>, public MxCrsMatrix<Scalar> {
  public:
    explicit MxYeeFitMu(RCP<MxEMSim<DIM> > theSim, bool invert = false);
    
    virtual const MxGridField<DIM> * getDomainField() const {return hfield;}

    virtual const MxGridField<DIM> * getRangeField() const {return hfield;}

    virtual const MxGrid<DIM> * getGrid() const {return grid;}

    void setTupleCoeffs(std::vector<double> coeffs);

    double * getTupleCoeffsView();

    int numTupleCoeffs() const {return numCoeffs;}

    std::vector<double> getTupleCoeffsCopy() const;

    //virtual void setMatrix();

    //virtual const Epetra_CrsMatrix * getMatrixPtr() const {return matrix.get();}

    //virtual Epetra_CrsMatrix getMatrixCopy() const {return *matrix;}

    //virtual void clearMatrix() {matrix = Teuchos::null;}

    RCP<MxCrsMatrix<Scalar> > cellAveInvMuOperator() const;

    /**
     *   | / | /
     *   |/__|/
     *  /|  /|
     * / | / |
     *
     *  Orderings:
     *    [E_xij, E_yij-1k, E_yijk, E_yi+1j-1k, E_yi+1jk,
     *            E_zijk-1, E_zijk, E_yi+1jk-1, E_zi+1jk]
     */
    class Stencil {
      public:
        std::vector<Scalar> values;

        std::vector<MxIndex> columns;

        std::vector<MxComplex> bcfactors;

        std::vector<MxDimVector<int, DIM> > cells;

        std::vector<size_t> comps;
    };

    class Tuple {
      public:
        MxDimVector<MxIndex, DIM> columns;

        MxDimVector<MxDimVector<int, DIM>, DIM> cells;

        MxDimVector<size_t, DIM> comps;

        MxDimVector<size_t, DIM> stencilInds;
    };

    Stencil getStencil(size_t comp0, MxDimVector<int, DIM> cell) const;

    void getTuples(Stencil const & stencil, std::vector<Tuple> & tuples) const;

    MxDimVector<double, DIM> getNormal(Stencil const & stencil) const;

    MxDimMatrix<MxComplex, DIM> tupleUpdate(Tuple const & tuple,
      MxDimVector<double, DIM> normal);

  private:

    bool mInvert;

    RCP<MxEMSim<DIM> > sim;

    const MxGridField<DIM> * hfield;

    const MxGridField<DIM> * bfield;

    const MxGrid<DIM> * grid;

    bool epsIsDiagonal(const MxDimMatrix<double, 3> & epsTensor) const;

    void setMatrix();

    void setMatrix2dTE();

    int numCoeffs;

};


#endif
