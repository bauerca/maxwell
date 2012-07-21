#ifndef MX_CONVEX_POLYGON
#define MX_CONVEX_POLYGON

#include <vector>

#include "MxDimVector.hpp"
#include "MxDynMatrix.hpp"
#include "MxPolytope.hpp"


template<size_t DIM>
class MxConvexPolygon : public MxPolytope<DIM> {
  public:

    virtual double volumeFraction(const MxShape<DIM> & aShape, const MxDimVector<double, DIM> & midPt) const;

    virtual double volumeFraction(const MxShape<DIM> & aShape, 
const std::vector<MxDimVector<double, DIM> > & verts,
const std::vector<double> & fVals,
const std::vector<size_t> & edgeXInds,
const std::vector<MxDimVector<double, DIM> *> & edgeX) const;

    virtual double volume() const = 0;

    //virtual MxPolytope<DIM> intersection(const MxShape<DIM> & aShape) const;

    virtual std::vector<MxDimVector<double, DIM> > getVertices(const MxDimVector<double, DIM> & p) const = 0;

    virtual MxDimVector<double, DIM> getCorner(const MxShape<DIM> & aShape, 
const std::vector<MxDimVector<double, DIM> > & verts,
const std::vector<size_t> & edgeXInds,
const std::vector<MxDimVector<double, DIM> *> & edgeX) const;

    virtual MxDynMatrix<size_t> getEdgeVertMatrix() const {return edgeVertMatrix;}

    virtual ~MxConvexPolygon() {}


  protected:

    MxDynMatrix<size_t> edgeVertMatrix;

    virtual void setLists() {this->setEdgeVertsLists();}
    virtual void setEdgeVertsLists();


};

#endif
