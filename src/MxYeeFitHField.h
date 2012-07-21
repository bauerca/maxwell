#ifndef MX_YEE_FIT_H_FIELD
#define MX_YEE_FIT_H_FIELD

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
#include "MxYeeMagFieldBase.h"

#include "Teuchos_ParameterList.hpp"

template<size_t> class MxGridFieldIter;

class MxMap;

template<size_t DIM>
class MxYeeFitHField : public MxYeeMagFieldBase<DIM> {
  public:
    MxYeeFitHField(const MxGrid<DIM> * aGrid, MxPolType polarization);

    virtual const MxPolytope<DIM> & getCompPolytope(size_t comp) const {return *this->compPtopes[comp];}

#if 0
    void setBCs(Teuchos::ParameterList theBCList);

    virtual bool useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const;

    virtual size_t globCompIndx(size_t comp, const MxDimVector<int, DIM> & cell) const;

    virtual MxComplex getCompFactor(size_t comp, MxDimVector<int, DIM> cell) const;

    virtual double calcCompFrac(size_t comp, MxDimVector<int, DIM> cell, const MxShape<DIM> & aShape) const;

    MxDimVector<int, DIM> getInteriorComp(size_t comp, MxDimVector<int, DIM> cell) const;
#endif

    friend class MxGridFieldIter<DIM>;

#if 0
  private:
    std::string pol;

    MxDimVector<MxBCType, DIM> lBCs, uBCs;

    MxDimVector<double, DIM> phaseShifts;

    MxDimVector<bool, DIM> useLowerComps, useUpperComps;

    MxDimVector<int, DIM> gridRes;

    Teuchos::ParameterList bcList;

    // a component direction is always 3D, even in a sub-3D simulation
    std::vector<MxDimVector<double, 3> > vecCompDirs;

    std::vector<int> compDirs;
#endif
};


#endif

