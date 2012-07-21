#ifndef MX_SEGMENT
#define MX_SEGMENT

#include "MxDimVector.hpp"
#include "MxDynMatrix.hpp"
#include "MxPolytope.hpp"


template<size_t DIM>
class MxSegment : public MxPolytope<DIM> {
  public:
    MxSegment(const MxDimVector<double, DIM> & p1, const MxDimVector<double, DIM> & p2);

    MxSegment(double aLen, const MxDimVector<double, DIM> & aDir);

    ~MxSegment() {}

    virtual double volumeFraction(const MxShape<DIM> & aShape, const MxDimVector<double, DIM> & midPt) const;

    virtual double volumeFraction(const MxShape<DIM> & aShape, 
const std::vector<MxDimVector<double, DIM> > & verts,
const std::vector<double> & fVals,
const std::vector<size_t> & edgeXInds,
const std::vector<MxDimVector<double, DIM> *> & edgeX) const;

    virtual double volume() const {return len;}

    //virtual MxPolytope<DIM> intersection(const MxShape<DIM> & aShape) const;
    
    virtual std::vector<MxDimVector<double, DIM> > getVertices(const MxDimVector<double, DIM> & midPt) const;

    virtual void getEdgeX(const MxShape<DIM> & aShape,
const std::vector<MxDimVector<double, DIM> > & verts,
const std::vector<double> & fVals,
std::vector<size_t> & edgeXInds,
std::vector<MxDimVector<double, DIM> *> & edgeX) const;

    virtual MxDynMatrix<size_t> getEdgeVertMatrix() const {return edgeVertMatrix;}



  protected:

    size_t numVerts;

    double len;

    MxDimVector<double, DIM> dir;

    MxDynMatrix<size_t> edgeVertMatrix;


};

#endif
