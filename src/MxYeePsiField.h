#ifndef MX_YEE_PSI_FIELD
#define MX_YEE_PSI_FIELD

#include <string>
#include <vector>
#include <utility>

#include "MxDimVector.hpp"
#include "MxShape.hpp"
#include "MxCartSeg.hpp"
#include "MxPoint.hpp"
#include "MxGridField.hpp"
#include "MxYeeFitBField.h"
#include "MxGrid.h"
#include "MxTypes.h"

#include "Teuchos_ParameterList.hpp"

template<size_t> class MxGridFieldIter;

class MxMap;

template<size_t DIM>
class MxYeePsiField : public MxGridField<DIM> {
  public:
    MxYeePsiField(const MxGrid<DIM> * aGrid, const MxYeeFitBField<DIM> * bfield);

    virtual const MxPolytope<DIM> & getCompPolytope(size_t comp) const {
      return *this->compPtopes[comp];
    }

    virtual const MxDimVector<double, DIM> & getCompCellCoord(size_t comp) const {return this->compCellCoords[comp];}

    void setBCs(MxDimVector<MxBCType, DIM> lower,
      MxDimVector<MxBCType, DIM> upper);

    virtual bool useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const;

    virtual MxComplex getCompFactor(size_t comp, MxDimVector<int, DIM> cell) const;

#if 0
    void setBCs(Teuchos::ParameterList theBCList);

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

    MxDimVector<int, DIM> gridRes;

    Teuchos::ParameterList bcList;
#endif

};


#endif

