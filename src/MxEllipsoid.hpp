#ifndef MX_ELLIPSOID
#define MX_ELLIPSOID

// Ellipsoid is a 3D object. If you want an ellipse, use Ellipse!

#include <string>

#include "MxShape.hpp"
#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxUtil.hpp"
#include "MxSphere.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLObject.hpp"

template<size_t DIM>
class MxEllipsoid : public MxShape<DIM> {
  public:
    MxEllipsoid(const MxDimVector<double, DIM> & location, const MxDimVector<double, 3> & axes, bool insideOut = false);

    // ctor from XML node
    MxEllipsoid(const Teuchos::XMLObject & node);

    virtual MxEllipsoid<DIM> * clone() const {return new MxEllipsoid(*this);}

    virtual bool hasGradFunc() const {return true;}

    virtual typename MxShape<DIM>::DimVecDblPair boundingBox() const {return bbox;}

  protected:
    
    double _func(const MxVecD3 & p) const {return 1.0 - p.dot(invAxes2 * p);}

    MxVecD3 _gradFunc(const MxVecD3 & p) const {return -2.0 * invAxes2 * p;}

  private:

    MxVecD3 invAxes2;

    typename MxShape<DIM>::DimVecDblPair bbox;
};


template<size_t DIM>
MxEllipsoid<DIM>::MxEllipsoid(const MxDimVector<double, DIM> & location,
const MxVecD3 & axes, bool insideOut) : 
invAxes2(1.0 / (axes * axes)) {
  MxShape<DIM>::init("ellipsoid");
  this->translate(MxVecD3(location));
}

template<size_t DIM>
MxEllipsoid<DIM>::MxEllipsoid(const Teuchos::XMLObject & node) {
  MxShape<DIM>::init(node);

  MxVecD3 axes;
  axes.strFill(MxUtil::XML::getAttr("axes", node));
  invAxes2 = 1.0 / (axes * axes);
}

/*
// gross over-estimation of the bounding box. This is just a cube based on
// the longest ellipsoidal axis.
template<size_t DIM>
void MxEllipsoid<DIM>::setBBox() {
  double maxAxis = MxDimVector<double, DIM>(axes).max();
  bbox.first = location - MxDimVector<double, 3>(maxAxis);
  bbox.second = location + MxDimVector<double, 3>(maxAxis);
}

// Defining function
template<size_t DIM>
inline double MxEllipsoid<DIM>::func(const MxDimVector<double, DIM> & p) const {

  // passive rotation of coordinates, actively rotates shape:
  MxDimVector<double, 3> pr(p - location);
  pr = RT * pr;

  pr /= axes; // scaling by ellipsoid axes
  return 1. - pr.dot(pr);
}

// Gradient of the boundary function for getting the normal
// to the surface.
template<size_t DIM>
inline MxDimVector<double, DIM> MxEllipsoid<DIM>::gradFunc (const MxDimVector<double, DIM> & p) const {
  // passive rotation of coordinates, actively rotates shape:
  MxDimVector<double, 3> pr(p - location);
  pr = RT * pr;
  pr /= -0.5 * axes * axes;
  return MxDimVector<double, DIM>(R * pr);
}
*/

#endif
