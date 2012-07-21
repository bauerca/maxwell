#ifndef MX_CONVEX_POLYHEDRON
#define MX_CONVEX_POLYHEDRON

#include <cmath>

#include "MxDimVector.hpp"
#include "MxDynMatrix.hpp"
#include "MxPolytope.hpp"


class MxConvexPolyhedron : public MxPolytope<3> {
  public:

    virtual double volumeFraction(const MxShape<3> & aShape, const MxDimVector<double, 3> & midPt) const;

    virtual double volumeFraction(const MxShape<3> & aShape, 
const std::vector<MxDimVector<double, 3> > & verts,
const std::vector<double> & fVals,
const std::vector<size_t> & edgeXInds,
const std::vector<MxDimVector<double, 3> *> & edgeX) const;

    // this should be simple enough for derived classes to implement
    virtual double volume() const = 0;

    //virtual MxPolytope<DIM> intersection(const MxShape<DIM> & aShape) const;
    
    // still purely virtual!
    virtual std::vector<MxDimVector<double, 3> > getVertices(const MxDimVector<double, 3> & p) const = 0;

    virtual MxDynMatrix<size_t> getEdgeVertMatrix() const {return edgeVertMatrix;}

    virtual MxDynMatrix<size_t> getEdgeFaceMatrix() const {return edgeFaceMatrix;}

    virtual ~MxConvexPolyhedron() {}


  protected:

    size_t numFaces;

    MxDynMatrix<size_t> edgeVertMatrix;
    MxDynMatrix<size_t> edgeFaceMatrix;
    MxDynMatrix<size_t> faceVertMatrix;

    std::vector<std::vector<size_t> > edgeFacesLists;
    std::vector<std::vector<size_t> > faceEdgesLists;
    std::vector<std::vector<size_t> > faceVertsLists;

    virtual void setEdgeVertsLists();
    virtual void setEdgeFacesLists();
    virtual void setFaceEdgesLists();
    virtual void setFaceVertsLists();

    virtual void setLists();


};

#endif
