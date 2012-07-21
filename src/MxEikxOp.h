#ifndef MX_EIKX_OP
#define MX_EIKX_OP

#include "MxTypes.h"
#include "MxDimVector.hpp"
#include "MxGridFieldFuncOp.h"

template<size_t DIM>
class MxEikxOp : public MxGridFieldFuncOp<DIM, MxComplex> {
  public:
    MxEikxOp(MxDimVector<double, DIM> blochK, double sign, double shift,
        const MxGridField<DIM> * aGridField,
        const MxGrid<DIM> * aGrid);

    virtual ~MxEikxOp() {};

  protected:

    virtual MxComplex func(size_t comp, MxDimVector<int, DIM> cell) const;

  private:
    
    std::complex<double> I;

    double mSign;

    double mShift;

    MxDimVector<double, DIM> mBlochK;

};

#endif
