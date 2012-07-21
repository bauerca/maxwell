#ifndef MX_TORUS
#define MX_TORUS

// Torus is a 3D object. If you want an annulus, use MxAnnulus!

#include <string>

#include "MxUtil.hpp"
#include "MxShape.hpp"
#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"

#include "Teuchos_XMLObject.hpp"

template<size_t DIM>
class MxTorus : public MxShape<DIM> {
  public:
    // axis is always 3D.
    MxTorus(const MxDimVector<double, DIM> & midPt, const MxDimVector<double, 3> & axis, double majorRadius, double minorRadius);

    MxTorus(const Teuchos::XMLObject & node);

    virtual MxTorus<DIM> * clone() const {return new MxTorus<DIM>(*this);}

    virtual bool hasGradFunc() const {return true;}

    virtual typename MxShape<DIM>::DimVecDblPair boundingBox() const {return bbox;}

  protected:
    
    virtual double _func(const MxDimVector<double, 3> & p) const {
      return r2 - pow(R - (complAxisProj * p).norm(), 2) - pow(axis.dot(p), 2);
    }

    virtual MxDimVector<double, 3> _gradFunc(const MxDimVector<double, 3> & p) const {
      return 2.0 * ((1.0 - R / (complAxisProj * p).norm()) * complAxisProj * p - axis.dot(p) * axis);
    }

  private:

    MxDimVector<double, 3> axis;

    double R, r, r2;

    // complementary axis projection
    MxDimMatrix<double, 3> complAxisProj;

    typename MxShape<DIM>::DimVecDblPair bbox;
};

 
template<size_t DIM>
MxTorus<DIM>::MxTorus(const MXDIMVEC & midPt, const MxVecD3 & torusAxis,
double majorRadius, double minorRadius) : 
axis(torusAxis / torusAxis.norm()), R(majorRadius), r(minorRadius),
complAxisProj(MxUtil::eye<double, 3>()) {
  MxShape<DIM>::init("torus");

  r2 = r * r;

  complAxisProj -= MxDimMatrix<double, 3>(axis, axis);

  // simple bounding box (gross overestimation of true bounding box)
  bbox.first = midPt - MxDimVector<double, DIM>(R + r);
  bbox.second = midPt + MxDimVector<double, DIM>(R + r);
}

template<size_t DIM>
inline
MxTorus<DIM>::MxTorus(const Teuchos::XMLObject & node) :
complAxisProj(MxUtil::eye<double, 3>()) {
  MxShape<DIM>::init(node);  

  R = atof(MxUtil::XML::getAttr("major radius", node).c_str());
  r = atof(MxUtil::XML::getAttr("minor radius", node).c_str());
  r2 = r * r;

  axis.strFill(MxUtil::XML::getAttr("axis", node));
  axis /= axis.norm();
  complAxisProj -= MxDimMatrix<double, 3>(axis, axis);
}

#endif
