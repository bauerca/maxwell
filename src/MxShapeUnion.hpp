///////////////////////////
//
//  MxShapeUnion represents the union of all the subshapes given to
//  this class. The union shape function is formed by multiplying all the subshape
//  functions together, then manipulating the resulting sign (positive if point
//  is in any shape). There is a small problem: in Mx, 0 is considered OUTSIDE
//  a shape; since 0s can occur INSIDE the shape union (where the boundary of
//  one subshape lies inside another subshape), a point that is actually inside
//  the shape union can be considered by Mx to be outside. To fix this, when a point, x,
//  lies inside any subshape (that is, f_i(x) > 0 for subshape i), then all f_j(x) = 0
//  for j \neq i are ignored in the product of subshape function values.
//
//

#ifndef MX_SHAPE_UNION
#define MX_SHAPE_UNION

#include <limits>

#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxShape.hpp"

template<size_t DIM>
class MxShapeUnion : public MxShape<DIM> {
  public:
    MxShapeUnion() : numSubShapes(0) {
      MxShape<DIM>::init("union");
    }

    int add(const MxShape<DIM> * aShape) {
      shapes.push_back(aShape);
      numSubShapes += aShape->getNumSubShapes();
      return numSubShapes;
    }

    virtual MxShapeUnion<DIM> * clone() const {return new MxShapeUnion<DIM>(*this);}

    virtual bool hasGradFunc() const {return true;}

    virtual typename MxShape<DIM>::DimVecDblPair boundingBox() const {return bbox;}

    virtual int getNumSubShapes() const {return numSubShapes;}

    virtual MxShape<DIM> * getSubShape(int indx) const;

  protected:

    virtual double _func(const MxVecD3 & p) const;

    virtual double _func(const MxVecD3 & p, int & indx) const;

    virtual MxVecD3 _gradFunc(const MxVecD3 & p) const;

  private:
    int numSubShapes;

    std::vector<const MxShape<DIM> *> shapes;

    typename MxShape<DIM>::DimVecDblPair bbox;
};


// if inside, choose the largest, positive shape func value
// if outside, choose the smallest (in magnitude) shape func value
template<size_t DIM>
inline
double MxShapeUnion<DIM>::_func(const MxVecD3 & p) const {
  MxDimVector<double, DIM> pp(p);

  double f;
  double fMax = -std::numeric_limits<double>::max();

  typename std::vector<const MxShape<DIM> *>::const_iterator iter;

  for (iter = shapes.begin(); iter != shapes.end(); ++iter) {
    f = (*iter)->func(pp);
    if (f > fMax) fMax = f;
  }

  return fMax;
}

// also need to implement version that provides sub-shape index
template<size_t DIM>
inline
double MxShapeUnion<DIM>::_func(const MxVecD3 & p, int & indx) const {
  MxDimVector<double, DIM> pp(p);

  double f;
  double fMax = -std::numeric_limits<double>::max();
  int tmpIndx, offset = 0;

  typename std::vector<const MxShape<DIM> *>::const_iterator iter;

  //std::cout << "shape union checking " << shapes.size() << " shapes\n";
  for (iter = shapes.begin(); iter != shapes.end(); ++iter) {
    //std::cout << "  shape union iter " << iter - shapes.begin() << "\n";
    f = (*iter)->func(pp, tmpIndx);
    if (f > fMax) {
      fMax = f;
      indx = offset + tmpIndx;
    }
    offset += (*iter)->getNumSubShapes();
  }
  //std::cout << "shape union done.\n";

  return fMax;
}


template<size_t DIM>
inline
MxShape<DIM> * MxShapeUnion<DIM>::getSubShape(int indx) const {
  int lb = 0, ub = 0;
  MxShape<DIM> * sh = 0;
  typename std::vector<const MxShape<DIM> *>::const_iterator iter;
  for (iter = shapes.begin(); iter != shapes.end(); ++iter) {
    ub += (*iter)->getNumSubShapes();
    if ((indx >= lb) and (indx < ub)) {
      sh = (*iter)->getSubShape(indx - lb);
      MxShape<DIM>::applyTransforms(*sh);
      return sh;
    }
    else
      lb = ub;
  }
  std::cout << "MxShapeUnion:getSubShape(...): Invalid sub-shape index!!\n";
  throw 1;
}



template<size_t DIM>
inline
MxVecD3 MxShapeUnion<DIM>::_gradFunc(const MxVecD3 & p) const {
  MxDimVector<double, DIM> pp(p);

  double f;
  double fMax = -std::numeric_limits<double>::max();

  typename std::vector<const MxShape<DIM> *>::const_iterator iter;
  const MxShape<DIM> * gradShape = 0;

  for (iter = shapes.begin(); iter != shapes.end(); ++iter) {
    f = (*iter)->func(pp);
    if (f > fMax) {
      fMax = f;
      gradShape = *iter;
    }
  }

  return MxVecD3(gradShape->gradFunc(pp));
}

#endif
