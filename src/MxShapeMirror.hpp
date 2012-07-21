

#ifndef MX_SHAPE_MIRROR
#define MX_SHAPE_MIRROR

#include <limits>

#include "MxUtil.hpp"
#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxShape.hpp"
#include "MxHalfSpace.hpp"

//! \class MxShapeMirror
  
/*! 
 * This is a wrapper class for an MxShape that will return the original shape plus its
   mirror image across a plane defined by the user.
*/
template<size_t DIM>
class MxShapeMirror : public MxShape<DIM> {
  public:
    /*!
     * Constructor:
     * \param theShape - A pointer to the shape being mirrored. This class will not destroy the target shape.
     * \param normal - The direction normal to the plane across which the shape will be mirrored (can
     * be of non-unit length)
     * \param pointInPlane - A point in the mirror plane (defaults to the origin)
     */
    MxShapeMirror(const MxShape<DIM> * theShape, const MxDimVector<double, DIM> & normal, const MxDimVector<double, DIM> & pointInPlane = 0);

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

    virtual MxShapeMirror<DIM> * clone() const {return new MxShapeMirror(*this);}

    virtual int getNumSubShapes() const {return 2 * numSubShapes;}

    virtual MxShape<DIM> * getSubShape(int indx) const;

    virtual bool hasGradFunc() const {return true;}

    virtual typename MxShape<DIM>::DimVecDblPair boundingBox() const {return bbox;}

  protected:

    virtual double _func(const MxVecD3 & p) const;

    virtual double _func(const MxVecD3 & p, int & indx) const;

    virtual MxVecD3 _gradFunc(const MxVecD3 & p) const;

  private:

    const MxShape<DIM> * shape;

    int numSubShapes;

    // Must create mirrored part of subshapes in order to return them via getSubShape(...). 
    // We'll use ShapeUnions as containers for the mirrored subshapes, which are mirrored by
    // applying a reflection to the ShapeUnion.
    //std::vector<Teuchos::RCP<MxShape<DIM> > > mirrorSubShapes;
    Teuchos::RCP<MxShape<DIM> > mirrorShape;

    // The mirror plane
    Teuchos::RCP<MxHalfSpace<3> > plane;

    typename MxShape<DIM>::DimVecDblPair bbox;
};


template<size_t DIM>
inline
MxShapeMirror<DIM>::MxShapeMirror(const MxShape<DIM> * theShape,
const MxDimVector<double, DIM> & normal,
const MxDimVector<double, DIM> & pointInPlane) :
shape(theShape), numSubShapes(theShape->getNumSubShapes()) {
  MxShape<DIM>::init("mirror");

  MxDimVector<double, 3> n(normal / normal.norm()), p(pointInPlane);
  plane = Teuchos::rcp(new MxHalfSpace<3>(p, n));

  MxShape<DIM> * sh = shape->clone();
  sh->reflect(n, p);
  mirrorShape = Teuchos::rcp(sh);
}


template<size_t DIM>
inline
double MxShapeMirror<DIM>::_func(const MxVecD3 & p) const {

  // if positive, do nothing, if negative, p is on the dark side of the mirror so
  // we use the reflected point
  double f = plane->func(p);

  if (f > 0.0)
    return shape->func(MxDimVector<double, DIM>(p));
  else
    return mirrorShape->func(MxDimVector<double, DIM>(p));
}

template<size_t DIM>
inline
double MxShapeMirror<DIM>::_func(const MxVecD3 & p, int & indx) const {

  // if positive, do nothing, if negative, p is on the dark side of the mirror so
  // we use the reflected point
  double f = plane->func(p);

  if (f > 0.0) {
    return shape->func(MxDimVector<double, DIM>(p), indx);
  }
  else {
    f = mirrorShape->func(MxDimVector<double, DIM>(p), indx);
    indx += numSubShapes;
    return f;
  }
}

template<size_t DIM>
inline
MxShape<DIM> * MxShapeMirror<DIM>::getSubShape(int indx) const {
  MxShape<DIM> * sh;
  if ((indx >= 0) and (indx < numSubShapes)) {
    sh = shape->getSubShape(indx);
    MxShape<DIM>::applyTransforms(*sh);
    return sh;
  }
  else if ((indx >= numSubShapes) and (indx < 2 * numSubShapes)) {
    sh = mirrorShape->getSubShape(indx - numSubShapes);
    MxShape<DIM>::applyTransforms(*sh);
    return sh;
  }
  else {
    std::cout << "MxShapeMirror::getSubShape(...): Invalid sub-shape index!!\n";
    throw 1;
  }
}

template<size_t DIM>
inline
MxVecD3 MxShapeMirror<DIM>::_gradFunc(const MxVecD3 & p) const {
  // if positive, do nothing, if negative, p is on the dark side of the mirror so
  // we use the reflected point
  double f = plane->func(p);

  if (f > 0.0)
    return MxVecD3(shape->gradFunc(MxDimVector<double, DIM>(p)));
  else
    return MxVecD3(mirrorShape->gradFunc(MxDimVector<double, DIM>(p)));
}

#endif
