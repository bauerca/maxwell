#ifndef MX_SPHERE
#define MX_SPHERE

#include <string>

#include "MxShape.hpp"
#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxUtil.hpp"

#include "Teuchos_XMLObject.hpp"

template<size_t DIM>
class MxSphere : public MxShape<DIM> {
  public:
    MxSphere(const MxDimVector<double, DIM> & location, double radius);

    // ctor from XML node
    MxSphere(const Teuchos::XMLObject & node);

    virtual MxSphere<DIM> * clone() const {return new MxSphere<DIM>(*this);}
    
    void setRotation (const MxDimVector<double, 3> & angles, const std::string & order);

    virtual bool hasGradFunc() const {return true;}

    virtual typename MxShape<DIM>::DimVecDblPair boundingBox() const {return bbox;}

  protected:

    virtual double _func(const MxDimVector<double, 3> & p) const {return 1.0 - p.dot(p) / r2;}

    virtual MxDimVector<double, 3> _gradFunc(const MxDimVector<double, 3> & p) const {return -2.0 * p / r2;}

  private:

    double r;

    double r2;

    typename MxShape<DIM>::DimVecDblPair bbox;
};
 
template<size_t DIM>
MxSphere<DIM>::MxSphere(const MxDimVector<double, DIM> & location, double radius) :
r(radius), r2(radius * radius) {
  MxShape<DIM>::init("sphere");
  this->translate(MxVecD3(location));
}


template<size_t DIM>
MxSphere<DIM>::MxSphere(const Teuchos::XMLObject & node) {
  MxShape<DIM>::init(node);
  if (node.getTag() == "Sphere") {
    r = atof(MxUtil::XML::getAttr("radius", node, "1").c_str());
    r2 = r * r;
  }
  else {
    std::cout << "Tag is not 'Sphere'!";
    throw 1;
  }
}


#endif
