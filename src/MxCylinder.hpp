#ifndef MX_CYLINDER
#define MX_CYLINDER

// Cylinder is a 3D object.

#include <string>

#include "MxUtil.hpp"
#include "MxShape.hpp"
#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"

#include "Teuchos_XMLObject.hpp"

template<size_t DIM>
class MxCylinder : public MxShape<DIM> {
  public:
    MxCylinder(const MxDimVector<double, DIM> & p, const MxDimVector<double, 3> & axis, double radius);

    MxCylinder(const Teuchos::XMLObject & node);

    virtual MxCylinder<DIM> * clone() const;

    virtual bool hasGradFunc() const {return true;}

    virtual typename MxShape<DIM>::DimVecDblPair boundingBox() const {return bbox;}

#if 0
    bool boundaryPoint(const MxDimVector<double, DIM> & p1, const MxDimVector<double, DIM> & p2, MxDimVector<double, DIM> & bndryPt) const;
#endif

  protected:

    virtual double _func(const MxDimVector<double, 3> & p) const {
      return r2 - p.dot(complAxisProj * p);
    }

    virtual MxDimVector<double, 3> _gradFunc(const MxDimVector<double, 3> & p) const {
      return -2.0 * complAxisProj * p;
    }

  private:

    MxVecD3 axis;

    double r, r2;

    // complementary axis projection
    MxDimMatrix<double, 3> complAxisProj;

    typename MxShape<DIM>::DimVecDblPair bbox;
};

template<size_t DIM> 
inline
MxCylinder<DIM>::MxCylinder(const MxDimVector<double, DIM> & midPt, const MxVecD3 & cylAxis,
double radius) : 
axis(cylAxis / cylAxis.norm()), r(radius), r2(radius * radius), 
complAxisProj(MxUtil::eye<double, 3>()) {
  MxShape<DIM>::init("cylinder");

  complAxisProj -= MxDimMatrix<double, 3>(axis, axis);

  // simple bounding box (gross overestimation of true bounding box)
  bbox.first = midPt - MxDimVector<double, DIM>(r);
  bbox.second = midPt + MxDimVector<double, DIM>(r);

  this->translate(MxVecD3(midPt));
}


template<size_t DIM> 
inline
MxCylinder<DIM>::MxCylinder(const Teuchos::XMLObject & node) :
complAxisProj(MxUtil::eye<double, 3>()) {
  MxShape<DIM>::init(node);

  r = atof(MxUtil::XML::getAttr("radius", node).c_str());
  r2 = r * r;

  axis.strFill(MxUtil::XML::getAttr("axis", node));
  axis /= axis.norm();
  complAxisProj -= MxDimMatrix<double, 3>(axis, axis);
}

template<size_t DIM>
inline
MxCylinder<DIM> * MxCylinder<DIM>::clone() const {return new MxCylinder<DIM>(*this);}

#if 0
template<>
inline
bool MxCylinder<3>::boundaryPoint(const MxDimVector<double, 3> & p1, const MxDimVector<double, 3> & p2, MxDimVector<double, 3> & bndryPt) const {
  
  MxDimVector<double, 3> dir(p1 - p2);
  dir /= dir.norm();

  MxDimVector<double, 3> prPerp(p1 - p0);
  prPerp -= dir.dot(prPerp) * dir;
  
  MxDimVector<double, 3> projPrPerp(complAxisProj * prPerp);
  MxDimVector<double, 3> projDir(complAxisProj * dir);

  double c1, c2, c3, c4;
  c1 = projDir.dot(prPerp);
  c2 = projDir.dot(dir);
  c3 = projPrPerp.dot(prPerp);
  double parRootPos = (-c1 + sqrt(c1 * c1 - c2 * (r2 - c3 * c3))) / c2;
  double parRootNeg = (-c1 - sqrt(c1 * c1 - c2 * (r2 - c3 * c3))) / c2;

  // now test which root is inside the given bounds

  MxDimVector<double, 3> rootPos, rootNeg; 
  rootPos = parRootPos * dir + prPerp + p0;
  rootNeg = parRootNeg * dir + prPerp + p0;
  
  if ((p1 - rootPos).dot(p2 - rootPos) < 0) {
    bndryPt = rootPos;
    return true;
  }
  else {
    bndryPt = rootNeg;
    return true;
  }
}
#endif

#endif
