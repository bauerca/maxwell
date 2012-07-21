#ifndef MX_POINT
#define MX_POINT

#include "MxPolytope.hpp"

template<size_t DIM>
class MxPoint : public MxPolytope<DIM> {
  public:
    MxPoint() {}

    virtual double volumeFraction(const MxShape<DIM> & aShape, const MxDimVector<double, DIM> & p) const {return aShape.func(p) > 0 ? 1 : 0;}

    virtual double volumeFraction(const MxShape<DIM> & aShape, 
const std::vector<MxDimVector<double, DIM> > & verts,
const std::vector<double> & fVals,
const std::vector<size_t> & edgeXInds,
const std::vector<MxDimVector<double, DIM> *> & edgeX) const {return fVals[0] > 0 ? 1 : 0;}

    virtual double volume() const {return 1;}

    virtual std::vector<MxDimVector<double, DIM> > getVertices(const MxDimVector<double, DIM> & p) const {return std::vector<MxDimVector<double, DIM> >(1, p);}

    virtual void getEdgeX(const MxShape<DIM> & aShape,
const std::vector<MxDimVector<double, DIM> > & verts,
const std::vector<double> & fVals,
std::vector<size_t> & edgeXInds,
std::vector<MxDimVector<double, DIM> *> & edgeX) const {;}
};

#endif
