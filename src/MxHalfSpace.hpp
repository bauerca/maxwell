#ifndef MX_HALF_SPACE
#define MX_HALF_SPACE

// HalfSpace: shape function is positive on one side of the plane, negative on the other

#include <string>

#include "MxUtil.hpp"
#include "MxShape.hpp"
#include "MxDimVector.hpp"

template<size_t DIM>
class MxHalfSpace : public MxShape<DIM> {
  public:
    MxHalfSpace(const MxDimVector<double, DIM> & pointInPlane, const MxDimVector<double, DIM> & normal) : 
    n(normal / normal.norm()) {
      MxShape<DIM>::init("plane");
      this->translate(MxVecD3(pointInPlane));
    }

    MxHalfSpace(const Teuchos::XMLObject & node);

    virtual MxHalfSpace<DIM> * clone() const {return new MxHalfSpace<DIM>(*this);}

    virtual bool hasGradFunc() const {return true;}

    virtual typename MxShape<DIM>::DimVecDblPair boundingBox() const {return bbox;}

  protected:
    
    virtual double _func(const MxVecD3 & p) const {return n.dot(p);}

    virtual MxVecD3 _gradFunc(const MxVecD3 & p) const {return n;}

  private:

    MxVecD3 n;

    typename MxShape<DIM>::DimVecDblPair bbox;
};
 
template<size_t DIM>
inline
MxHalfSpace<DIM>::MxHalfSpace(const Teuchos::XMLObject & node) {
  MxShape<DIM>::init(node);

  n.strFill(MxUtil::XML::getAttr("normal", node));
  n /= n.norm();
}

#endif
