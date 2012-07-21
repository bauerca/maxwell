#ifndef MX_POLYTOPE
#define MX_POLYTOPE

// Polytope
// A polytope is a collection of connected vertices. More important, a polytope
// is something of which you can calculate the volume (or the volume of its 
// intersection with a shape)
//

#include <utility>
#include <vector>

#include "MxDimVector.hpp"
#include "MxDynMatrix.hpp"
#include "MxShape.hpp"



//#include "MxGeoObj.h"

//template<size_t DIM> class MxShape;
//template<typename T, size_t DIM> class MxDimVector;

enum MxPolytopeState {
  IN,
  OUT,
  CUT1,
  ON
};

template<size_t DIM>
class MxPolytope {
  public:

    virtual double volume() const = 0;

    // get volume fraction for polytope at position p, w/out changing internal
    // location of polytope
    virtual double volumeFraction(const MxShape<DIM> & aShape, const MxDimVector<double, DIM> & p) const = 0;

    virtual double volumeFraction(const MxShape<DIM> & aShape, 
const std::vector<MxDimVector<double, DIM> > & verts,
const std::vector<double> & fVals,
const std::vector<size_t> & edgeXInds,
const std::vector<MxDimVector<double, DIM> *> & edgeX) const = 0;
    //virtual MxPolyTope<DIM> intersection(const MxShape<DIM> & aShape) const = 0;

    // return the vertices' locations given some midpoint, centroid, etc.
    virtual std::vector<MxDimVector<double, DIM> > getVertices(const MxDimVector<double, DIM> &) const = 0;

    virtual void getEdgeX(const MxShape<DIM> & aShape,
const std::vector<MxDimVector<double, DIM> > & verts,
const std::vector<double> & fVals,
std::vector<size_t> & edgeXInds,
std::vector<MxDimVector<double, DIM> *> & edgeX) const;

    virtual MxPolytopeState getState(const MxShape<DIM> & aShape, const std::vector<MxDimVector<double, DIM> > & verts, std::vector<double> & fVals) const;

    //virtual MxDynMatrix<size_t> getEdgeVertMatrix() const = 0;

    virtual ~MxPolytope() {}

  protected:

    size_t numVerts;
    size_t numEdges;

    std::vector<std::vector<size_t> > edgeVertsLists;
};

#endif
