#ifndef MX_CONE
#define MX_CONE

// Cone is a 3D object. If you want an annulus, use MxAnnulus!

#include <string>
#include <cmath>

#include "MxUtil.hpp"
#include "MxShape.hpp"
#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"

#include "Teuchos_XMLObject.hpp"

template<size_t DIM>
class MxCone : public MxShape<DIM> {
  public:
    MxCone(const MxDimVector<double, 3> & vertex, const MxDimVector<double, 3> & axis, double angle);
    
    MxCone(const Teuchos::XMLObject & node);

    virtual MxCone<DIM> * clone() const {return new MxCone<DIM>(*this);}

    virtual bool hasGradFunc() const {return true;}

    virtual typename MxShape<DIM>::DimVecDblPair boundingBox() const {return bbox;}

  protected:

    virtual double _func(const MxDimVector<double, 3> & p) const {
      return pow(tanTheta * axis.dot(p), 2) - p.dot(complAxisProj * p);
    }

    virtual MxDimVector<double, 3> _gradFunc(const MxDimVector<double, 3> & p) const {
      return 2.0 * (tanTheta * tanTheta * axis.dot(p) * axis - complAxisProj * p);
    }

  private:

    MxDimVector<double, 3> axis;

    double tanTheta;

    // complementary axis projection
    MxDimMatrix<double, 3> complAxisProj;

    typename MxShape<DIM>::DimVecDblPair bbox;
};

template<size_t DIM> 
inline
MxCone<DIM>::MxCone(const MxVecD3 & vertex, const MxVecD3 & coneAxis, double coneAngle) : 
axis(coneAxis / coneAxis.norm()), tanTheta(tan(coneAngle)), complAxisProj(MxUtil::eye<double, 3>()) {
  MxShape<DIM>::init("cone");

  complAxisProj -= MxDimMatrix<double, 3>(axis, axis);

  this->translate(vertex);
}

template<size_t DIM> 
inline
MxCone<DIM>::MxCone(const Teuchos::XMLObject & node) : 
complAxisProj(MxUtil::eye<double, 3>()) {
  MxShape<DIM>::init(node);

  tanTheta = tan(atof(MxUtil::XML::getAttr("angle", node).c_str()));

  axis.strFill(MxUtil::XML::getAttr("axis", node));
  axis /= axis.norm();
  complAxisProj -= MxDimMatrix<double, 3>(axis, axis);
}

#endif
