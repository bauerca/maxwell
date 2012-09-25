#ifndef MX_YEE_FIT_E_FIELD
#define MX_YEE_FIT_E_FIELD

#include <string>
#include <vector>
#include <utility>

#include "MxDimVector.hpp"
#include "MxShape.hpp"
#include "MxCartSeg.hpp"
#include "MxPoint.hpp"
#include "MxGridField.hpp"
#include "MxGrid.h"
#include "MxTypes.h"
#include "MxYeeElecFieldBase.h"
#include "MxYeeFitBField.h"

#include "Teuchos_ParameterList.hpp"

template<size_t> class MxGridFieldIter;

class MxMap;

template<size_t DIM>
class MxYeeFitEField : public MxYeeElecFieldBase<DIM> {
  public:
    // needs bfield for Dey Mittra
    MxYeeFitEField(const MxGrid<DIM> * aGrid, const MxYeeFitBField<DIM> * bfield,
        MxPolType polarization = TE);

    virtual const MxPolytope<DIM> & getCompPolytope(size_t comp) const {
      return *this->compPtopes[comp];
    }

    virtual bool useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const;

    virtual MxComplex getCompFactor(size_t comp, MxDimVector<int, DIM> cell) const;

    //void setBCs(Teuchos::ParameterList theBCList);

// these have been moved to super class (12/30/2011)
#if 0
    virtual bool useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const;

    virtual size_t globCompIndx(size_t comp, const MxDimVector<int, DIM> & cell) const;

    virtual MxComplex getCompFactor(size_t comp, MxDimVector<int, DIM> cell) const;

    virtual double calcCompFrac(size_t comp, MxDimVector<int, DIM> cell, const MxShape<DIM> & aShape) const;

    MxDimVector<int, DIM> getInteriorComp(size_t comp, MxDimVector<int, DIM> cell) const;
#endif

    friend class MxGridFieldIter<DIM>;

  private:
    
    const MxYeeFitBField<DIM> * mBfield;

#if 0
  private:

    MxDimVector<MxBCType, DIM> lBCs, uBCs;

    MxDimVector<double, DIM> phaseShifts;

    MxDimVector<bool, DIM> useLowerComps, useUpperComps;

    MxDimVector<int, DIM> gridRes;

    Teuchos::ParameterList bcList;
#endif
};


#endif

