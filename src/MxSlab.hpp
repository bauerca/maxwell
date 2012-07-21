#ifndef MX_SLAB
#define MX_SLAB

// MxSlab: shape function is positive on one side of the plane, negative on the other

#include <string>

#include "MxShape.hpp"
#include "MxHalfSpace.hpp"
#include "MxDimVector.hpp"

template<size_t DIM>
class MxSlab : public MxShape<DIM> {
  public:
    MxSlab(const MxDimVector<double, DIM> & pointInMidPlane, const MxDimVector<double, DIM> & normal, double thickness);

    MxSlab(const Teuchos::XMLObject & node);

    MxSlab(const MxSlab<DIM> & theSlab);

    MxSlab<DIM> & operator=(const MxSlab<DIM> & theSlab);

    virtual ~MxSlab();

    virtual MxSlab<DIM> * clone() const {return new MxSlab<DIM>(*this);}

    virtual bool hasGradFunc() const {return true;}

    virtual typename MxShape<DIM>::DimVecDblPair boundingBox() const {return bbox;}

  protected:
    
    virtual double _func(const MxVecD3 & p) const {
      return slab.func(MxDimVector<double, DIM>(p));
    }

    virtual MxVecD3 _gradFunc(const MxVecD3 & p) const {
      return MxVecD3(slab.gradFunc(MxDimVector<double, DIM>(p)));
    }

  private:

    MxHalfSpace<DIM> * plane1, * plane2;

    MxShapeIntersection<DIM> slab;

    typename MxShape<DIM>::DimVecDblPair bbox;

    void build(double l, const MXDIMVEC & n);

};
 
template<size_t DIM>
MxSlab<DIM>::MxSlab(const MxDimVector<double, DIM> & pointInMidPlane, const MxDimVector<double, DIM> & normal,
double thickness) {
  MxShape<DIM>::init("slab");
  build(thickness, normal / normal.norm());
  this->translate(MxVecD3(pointInMidPlane));
}

template<size_t DIM>
MxSlab<DIM>::MxSlab(const Teuchos::XMLObject & node) {
  MxShape<DIM>::init(node);

  MxDimVector<double, DIM> n;
  double l;

  n.strFill(MxUtil::XML::getAttr("normal", node));
  n /= n.norm();

  l = atof(MxUtil::XML::getAttr("thickness", node).c_str());

  build(l, n);
}

template<size_t DIM>
MxSlab<DIM> & MxSlab<DIM>::operator=(const MxSlab<DIM> & theSlab) {
  MxShape<DIM>::operator=(theSlab);
  bbox = theSlab.bbox;
  plane1 = new MxHalfSpace<DIM>(*theSlab.plane1);
  plane2 = new MxHalfSpace<DIM>(*theSlab.plane2);
  slab.add(plane1);
  slab.add(plane2);
}

template<size_t DIM>
MxSlab<DIM>::MxSlab(const MxSlab<DIM> & theSlab) : MxShape<DIM>(theSlab), bbox(theSlab.bbox) {
  plane1 = new MxHalfSpace<DIM>(*theSlab.plane1);
  plane2 = new MxHalfSpace<DIM>(*theSlab.plane2);
  slab.add(plane1);
  slab.add(plane2);
}

template<size_t DIM>
void MxSlab<DIM>::build(double l, const MxDimVector<double, DIM> & n) {
  plane1 = new MxHalfSpace<DIM>(-0.5 * l * n, n);
  plane2 = new MxHalfSpace<DIM>( 0.5 * l * n, -n);
  slab.add(plane1);
  slab.add(plane2);
}


template<size_t DIM>
MxSlab<DIM>::~MxSlab() {
  delete plane1;
  delete plane2;
}


#endif
