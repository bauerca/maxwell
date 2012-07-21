

#ifndef MX_SHAPE_REPEAT
#define MX_SHAPE_REPEAT

#include <limits>

#include "MxUtil.hpp"
#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxShape.hpp"

#include "Teuchos_RCP.hpp"

//! \class MxShapeRepeat
  
/*! 
 * This is a wrapper class for an MxShape that will return the original shape plus its
   mirror image across a plane defined by the user.
*/
template<size_t DIM>
class MxShapeRepeat : public MxShape<DIM> {
  public:
    /*!
     * Constructor:
     * \param theShape - A pointer to the shape being mirrored. This class will not destroy the target shape.
     * \param normal - The direction normal to the plane across which the shape will be mirrored (can
     * be of non-unit length)
     * \param pointInPlane - A point in the mirror plane (defaults to the origin)
     */
    MxShapeRepeat(const MxShape<DIM> * theShape, const MxDimVector<double, DIM> & origin, const MxDimVector<double, DIM> & direction, double step, int numPos, int numNeg);

    /*! 
     * We need to override the MxShape implementation of func because we need to do the coordinate
     * mirroring \e before the general Galilean transformations of the coordinates.
     * \param p - Shape function evaluation point
     */
    // no we don't
    //double func(const MxDimVector<double, DIM> & p) const;

    /*! 
     * We need to override the MxShape implementation of gradFunc because we need to do the coordinate
     * mirroring \e before the general Galilean transformations of the coordinates.
     * \param p - Shape function evaluation point
     */
    // no we don't
    //MxDimVector<double, DIM> gradFunc(const MxDimVector<double, DIM> & p) const;

    virtual MxShapeRepeat<DIM> * clone() const {return new MxShapeRepeat<DIM>(*this);}

    virtual bool hasGradFunc() const {return true;}

    virtual typename MxShape<DIM>::DimVecDblPair boundingBox() const {return bbox;}

    virtual int getNumSubShapes() const {return numSubShapes;}

    virtual MxShape<DIM> * getSubShape(int indx) const;

  protected:

    virtual double _func(const MxVecD3 & p) const;

    virtual double _func(const MxVecD3 & p, int & indx) const;

    virtual MxVecD3 _gradFunc(const MxVecD3 & p) const;

  private:

    const MxShape<DIM> * shape;

    MxVecD3 o, dir;

    double s;

    int numSubShapes;

    double np, nn;

    std::vector<Teuchos::RCP<MxShape<DIM> > > subShapes;

    typename MxShape<DIM>::DimVecDblPair bbox;
};


template<size_t DIM>
inline
MxShapeRepeat<DIM>::MxShapeRepeat(const MxShape<DIM> * theShape,
const MxDimVector<double, DIM> & origin,
const MxDimVector<double, DIM> & direction,
double step,
int numPos,
int numNeg) :
shape(theShape), o(origin), dir(direction / direction.norm()), s(step),
numSubShapes(theShape->getNumSubShapes() * (1 + numPos + numNeg)),
np(numPos), nn(-numNeg) {
  MxShape<DIM>::init("repeat");

  subShapes.resize(1 + numPos + numNeg);
  MxShape<DIM> * sh;
  for (int i = -numNeg; i <= numPos; ++i) {
    sh = shape->clone();
    sh->translate(double(i) * s * dir);
    subShapes[i + numNeg] = Teuchos::rcp(sh);
  }

}

template<size_t DIM>
inline
double MxShapeRepeat<DIM>::_func(const MxVecD3 & p) const {
  MxVecD3 pp;
  double slabPt = (p - o).dot(dir) / s + 0.5;
  if (slabPt >= 0.0) 
    pp = p - s * std::min(floor(slabPt), np) * dir;
  else
    pp = p - s * std::max(floor(slabPt), nn) * dir;

  return shape->func(MxDimVector<double, DIM>(pp));
}

template<size_t DIM>
inline
double MxShapeRepeat<DIM>::_func(const MxVecD3 & p, int & indx) const {
  MxVecD3 pp;
  //std::cout << "func eval point: "; p.print();
  double slabPt = (p - o).dot(dir) / s + 0.5;
  double dSlab;
  if (slabPt >= 0.0)
    dSlab = std::min(floor(slabPt), np);
  else
    dSlab = std::max(floor(slabPt), nn);

  pp = p - s * dSlab * dir;
  double res = shape->func(MxDimVector<double, DIM>(pp), indx);
  //std::cout << "dSlab - nn = " << dSlab - nn << "\n";
  indx = int(dSlab - nn) * shape->getNumSubShapes() + indx;
  return res;
}

template<size_t DIM>
inline
MxVecD3 MxShapeRepeat<DIM>::_gradFunc(const MxVecD3 & p) const {
  MxVecD3 pp;
  double slabPt = (p - o).dot(dir) / s + 0.5;
  if (slabPt >= 0.0) 
    pp = p - s * std::min(floor(slabPt), np) * dir;
  else
    pp = p - s * std::max(floor(slabPt), nn) * dir;

  return MxVecD3(shape->gradFunc(MxDimVector<double, DIM>(pp)));
}

template<size_t DIM>
inline
MxShape<DIM> * MxShapeRepeat<DIM>::getSubShape(int indx) const {
  int copy = indx / shape->getNumSubShapes();
  int copySubShapeIndx = indx % shape->getNumSubShapes();

  MxShape<DIM> * sh = subShapes[copy]->getSubShape(copySubShapeIndx);
  MxShape<DIM>::applyTransforms(*sh);
  return sh;
}


#endif
