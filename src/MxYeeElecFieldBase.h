#ifndef MX_YEE_ELEC_FIELD_BASE
#define MX_YEE_ELEC_FIELD_BASE

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

#include "Teuchos_ParameterList.hpp"

template<size_t> class MxGridFieldIter;

class MxMap;

/**
The YeeElecFieldBase class sets some GridField properties that are
common to all Yee electric fields. So far, these are just: 

  1) components are located on cell edges
  2) electromagnetic boundary conditions (PEC/PMC)
  3) specification of polarization for 2D simulations

The polytopes that are associated with each field component are
omitted in this class. Derived classes must set them.
*/
template<size_t DIM>
class MxYeeElecFieldBase : public MxGridField<DIM> {
  public:
    void initYeeElecFieldBase(const MxGrid<DIM> * aGrid,
      MxPolType polarization = TE);

    virtual const MxDimVector<double, DIM> & getCompCellCoord(size_t comp) const {return this->compCellCoords[comp];}

    //void setBCs(Teuchos::ParameterList theBCList);

    virtual void setBCs(MxDimVector<MxBCType, DIM> lower,
        MxDimVector<MxBCType, DIM> upper);

    virtual MxPolType getPolarization() const {return pol;}
// these have been moved to super class (12/30/2011)
#if 0
    virtual bool useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const;

    virtual size_t globCompIndx(size_t comp, const MxDimVector<int, DIM> & cell) const;

    virtual MxComplex getCompFactor(size_t comp, MxDimVector<int, DIM> cell) const;

    virtual double calcCompFrac(size_t comp, MxDimVector<int, DIM> cell, const MxShape<DIM> & aShape) const;

    MxDimVector<int, DIM> getInteriorComp(size_t comp, MxDimVector<int, DIM> cell) const;
#endif

    friend class MxGridFieldIter<DIM>;

  protected:
    MxPolType pol;

#if 0
    MxDimVector<MxBCType, DIM> lBCs, uBCs;

    MxDimVector<double, DIM> phaseShifts;

    MxDimVector<bool, DIM> useLowerComps, useUpperComps;

    MxDimVector<int, DIM> gridRes;

    Teuchos::ParameterList bcList;
#endif

    // a component direction is always 3D, even in a sub-3D simulation
    std::vector<MxDimVector<double, 3> > vecCompDirs;

    std::vector<int> compDirs;
};


#endif

