
#ifndef MX_SHAPE_SUBTRACT
#define MX_SHAPE_SUBTRACT

#include <limits>

#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxShape.hpp"
#include "MxShapeUnion.hpp"

template<size_t DIM>
class MxShapeSubtract : public MxShape<DIM> {
  public:
    MxShapeSubtract() : numSubShapes(0), baseShape(0) {
      MxShape<DIM>::init("subtract");
      rmShape.invert();
    }

    virtual MxShapeSubtract<DIM> * clone() const {return new MxShapeSubtract<DIM>(*this);}

    void setBaseShape(const MxShape<DIM> * aShape);

    void subtractShape(const MxShape<DIM> * aShape);

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

    // never owned by this object
    const MxShape<DIM> * baseShape;

    MxShapeUnion<DIM> rmShape;

    typename MxShape<DIM>::DimVecDblPair bbox;

};


template<size_t DIM>
inline
void MxShapeSubtract<DIM>::setBaseShape(const MxShape<DIM> * aShape) {
  baseShape = aShape;
  numSubShapes += aShape->getNumSubShapes();
  bbox = baseShape->boundingBox();
}

template<size_t DIM>
inline
void MxShapeSubtract<DIM>::subtractShape(const MxShape<DIM> * aShape) {
  rmShape.add(aShape);
  numSubShapes += aShape->getNumSubShapes();
}


template<size_t DIM>
inline
double MxShapeSubtract<DIM>::_func(const MxVecD3 & p) const {
  MxDimVector<double, DIM> pp(p);

  double fBase = baseShape->func(pp);
  bool inBase = fBase > 0;

  double fRm = rmShape.func(pp);
  bool inRm = fRm < 0;

  if (inRm and inBase)
    return fRm;
  else if (!inRm and !inBase)
    return fBase;
  else 
    return fBase < fRm ? fBase : fRm;
    // both are the same sign. If negative, we are inside the 
    // removal shape and outside the base shape (choose most negative). 
    // If positive, we are inside the base shape and outside the
    // removal shape (choose least positive).

}

template<size_t DIM>
inline
double MxShapeSubtract<DIM>::_func(const MxVecD3 & p, int & indx) const {
  MxDimVector<double, DIM> pp(p);

  int baseIndx;
  double fBase = baseShape->func(pp, baseIndx);
  bool inBase = fBase > 0;

  int rmIndx;
  double fRm = rmShape.func(pp, rmIndx);
  rmIndx += baseShape->getNumSubShapes();
  bool inRm = fRm < 0;

  if (inRm and inBase) {
    indx = rmIndx;
    return fRm;
  }
  else if (!inRm and !inBase) {
    indx = baseIndx;
    return fBase;
  }
  else if (fBase < fRm) {
    indx = baseIndx;
    return fBase;
  }
  else {
    indx = rmIndx;
    return fRm;
  }
    // both are the same sign. If negative, we are inside the 
    // removal shape and outside the base shape (choose most negative). 
    // If positive, we are inside the base shape and outside the
    // removal shape (choose least positive).
}


template<size_t DIM>
inline
MxShape<DIM> * MxShapeSubtract<DIM>::getSubShape(int indx) const {
  int lb = 0, ub = baseShape->getNumSubShapes();
  MxShape<DIM> * sh = 0;
  if ((indx >= lb) and (indx < ub))
    sh = baseShape->getSubShape(indx);

  lb = ub;
  ub += rmShape.getNumSubShapes();
  if ((indx >= lb) and (indx < ub))
    sh = rmShape.getSubShape(indx - lb);

  if (sh != 0) {
    MxShape<DIM>::applyTransforms(*sh);
    return sh;
  }

  std::cout << "MxShapeSubtract:getSubShape(...): Invalid sub-shape index!!\n";
  throw 1;
}


template<size_t DIM>
inline
MxVecD3 MxShapeSubtract<DIM>::_gradFunc(const MxVecD3 & p) const {
  MxDimVector<double, DIM> pp(p);

  double fBase = baseShape->func(pp);
  bool inBase = fBase > 0;

  double fRm = rmShape.func(pp);
  bool inRm = fRm < 0;

  if (inRm and inBase)
    return MxVecD3(rmShape.gradFunc(pp));
  else if (!inRm and !inBase)
    return MxVecD3(baseShape->gradFunc(pp));
  else 
    return fBase < fRm ? MxVecD3(baseShape->gradFunc(pp)) : MxVecD3(rmShape.gradFunc(pp));
    // both are the same sign. If negative, we are inside the 
    // removal shape and outside the base shape (choose most negative). 
    // If positive, we are inside the base shape and outside the
    // removal shape (choose least positive).
}

#endif
