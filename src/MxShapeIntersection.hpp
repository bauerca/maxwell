///////////////////////////
//
//  This class represents a pec. The shapes given to this
//  class are always combined by union.

#ifndef MX_SHAPE_INTERSECTION
#define MX_SHAPE_INTERSECTION

#include <utility>

#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxShape.hpp"

template<size_t DIM>
class MxShapeIntersection : public MxShape<DIM> {
  public:
    MxShapeIntersection() : numSubShapes(0) {
      MxShape<DIM>::init("intersection");
    }

    virtual MxShapeIntersection<DIM> * clone() const {return new MxShapeIntersection<DIM>(*this);}

    int add(const MxShape<DIM> * aShape);

    virtual bool hasGradFunc() const {return true;}

    virtual int getNumSubShapes() const {return numSubShapes;}

    virtual MxShape<DIM> * getSubShape(int indx) const;

    virtual typename MxShape<DIM>::DimVecDblPair boundingBox() const {return bbox;}

  protected:

    virtual double _func(const MxVecD3 & p) const;

    virtual double _func(const MxVecD3 & p, int & indx) const;

    virtual MxVecD3 _gradFunc(const MxVecD3 & p) const;

  private:
    int numSubShapes;

    std::vector<const MxShape<DIM> *> shapes;

    typename MxShape<DIM>::DimVecDblPair bbox;

    void setBBox();
};


template<size_t DIM>
inline
int MxShapeIntersection<DIM>::add(const MxShape<DIM> * aShape) {
  shapes.push_back(aShape);
  numSubShapes += aShape->getNumSubShapes();

  setBBox();

  return numSubShapes;
}

template<size_t DIM>
inline
void MxShapeIntersection<DIM>::setBBox() {
  // warning: if the added shape does not intersect any preexisting shapes,
  // the bounding box will be wrong. But really this situation should never
  // be encountered.

  if (numSubShapes == 1)
    bbox = shapes[0]->boundingBox();
  else {
    // otherwise the bounding box is the intersection of the underlying
    // shapes' bounding boxes

    bool lb1iInsideAll, ub1iInsideAll;
    int lb1i, ub1i, lb2i, ub2i;
    typename std::vector<const MxShape<DIM> *>::const_iterator iter1, iter2;

    for (iter1 = shapes.begin(); iter1 != shapes.end(); ++iter1) {
      // loop through each dimension of the simulation
      for (size_t i = 0; i < DIM; i++) {
        // get bounds for current shape in this dimension
        lb1i = (*iter1)->boundingBox().first[i];
        ub1i = (*iter1)->boundingBox().second[i];

        lb1iInsideAll = true;
        ub1iInsideAll = true;

        // loop through all other shapes. If one of the bounds for
        // the current shape lies within the bounds of all other shapes,
        // use it as a bound for the shape intersection
        for (iter2 = shapes.begin(); iter2 != shapes.end(); ++iter2) {
          // ignore self comparison
          if (iter1 == iter2)
            continue;
          else {
            lb2i = (*iter2)->boundingBox().first[i];
            ub2i = (*iter2)->boundingBox().second[i];

            // is current lower bound out-of-bounds of comparison shape?
            if (lb1i < lb2i or lb1i > ub2i)
              lb1iInsideAll = false;
            // is current upper bound out-of-bounds of comparison shape?
            if (ub1i < lb2i or ub1i > ub2i)
              ub1iInsideAll = false;
          }
        } // done looping over all other shapes

        if (lb1iInsideAll) bbox.first[i] = lb1i;
        if (ub1iInsideAll) bbox.second[i] = ub1i;
      }
    }
  }

}


template<size_t DIM>
inline
double MxShapeIntersection<DIM>::_func(const MxVecD3 & p) const {
  MxDimVector<double, DIM> pp(p);

  double f;
  double fMin = std::numeric_limits<double>::max();

  typename std::vector<const MxShape<DIM> *>::const_iterator iter;

  for (iter = shapes.begin(); iter != shapes.end(); ++iter) {
    f = (*iter)->func(pp);
    if (f < fMin) fMin = f;
  }

  return fMin;
}

template<size_t DIM>
inline
double MxShapeIntersection<DIM>::_func(const MxVecD3 & p, int & indx) const {
  MxDimVector<double, DIM> pp(p);

  double f;
  double fMin = std::numeric_limits<double>::max();
  int tmpIndx, offset = 0;

  typename std::vector<const MxShape<DIM> *>::const_iterator iter;

  //std::cout << "MxShapeIntersection::_func(..,indx)...\n";
  //std::cout << "shape intersection checking " << shapes.size() << " shapes\n";
  for (iter = shapes.begin(); iter != shapes.end(); ++iter) {
    //std::cout << "  shape intersection iter " << iter - shapes.begin() << "\n";
    f = (*iter)->func(pp, tmpIndx);
    //std::cout << "  " << (*iter)->getName() << " f = " << f << "\n";
    if (f < fMin) {
      fMin = f;
      indx = offset + tmpIndx;
    }
    offset += (*iter)->getNumSubShapes();
  }
  //std::cout << "shape intersection done.\n";

  return fMin;
}

template<size_t DIM>
inline
MxShape<DIM> * MxShapeIntersection<DIM>::getSubShape(int indx) const {
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
  std::cout << "MxShapeIntersection:getSubShape(...): Invalid sub-shape index!!\n";
  throw 1;
}

template<size_t DIM>
inline
MxVecD3 MxShapeIntersection<DIM>::_gradFunc(const MxVecD3 & p) const {
  MxDimVector<double, DIM> pp(p);

  double f;
  double fMin = std::numeric_limits<double>::max();

  typename std::vector<const MxShape<DIM> *>::const_iterator iter;
  const MxShape<DIM> * gradShape = 0;

  for (iter = shapes.begin(); iter != shapes.end(); ++iter) {
    f = (*iter)->func(pp);
    if (f < fMin) {
      fMin = f;
      gradShape = *iter;
    }
  }

  return MxVecD3(gradShape->gradFunc(pp));
}

#endif
