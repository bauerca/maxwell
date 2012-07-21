
#ifndef MX_SHAPE_DIFFERENCE
#define MX_SHAPE_DIFFERENCE

#include <vector>

#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxShape.hpp"

template<size_t DIM>
class MxShapeDifference : public MxShape<DIM> {
  public:
    MxShapeDifference() : numShapes(0) {}

    int add(const MxShape<DIM> * aShape) {
      shapes.push_back(aShape);
      numShapes++;
      return numShapes;
    }

    virtual double func(const MxDimVector<double, DIM> & p) const;

    virtual MxDimVector<double, DIM> gradFunc(const MxDimVector<double, DIM> & p) const;

    virtual bool hasGradFunc() const {return true;}

    int getNumShapes() const {return numShapes;}


  private:
    int numShapes;

    std::vector<const MxShape<DIM> *> shapes;

};

template<size_t DIM>
inline
double MxShapeDifference<DIM>::func(const MxDimVector<double, DIM> & p) const {
  double res = 1;
  double f;
  int numInside = 0;

  std::vector<const MxShape<DIM> *>::const_iterator iter;

  for (iter = shapes.begin(); iter != shapes.end(); ++iter) {
    f = (*iter)->func(p);
    if (f > 0) numInside++;
    if (f > fMax) fMax = f;
    res *= f;
  }

  return numInside == 1 ? fabs(res) : -fabs(res);
}

template<size_t DIM>
inline
MxDimVector<double, DIM> MxShapeDifference<DIM>::gradFunc(const MxDimVector<double, DIM> & p) const {

  int numInside = 0;
  double f, fprod = 1;

  // first get/save function values for all underlying shapes and calculate
  // the product
  double fVals[numShapes];
  for (int i = 0; i < numShapes; i++) {
    f = shapes[i]->func(p);
    fprod *= f;
    if (f > 0) numInside++;
    fVals[i] = f;
  }

  // when this->func(p) callously forces a sign-change on the product of
  // the underlying shape func values, the sign of the gradient must be
  // changed as well. If 'fprod' gives the same value as this->func(p), then
  // no sign change is required.
  double gradSign;
  if (fprod * (numInside == 1 ? fabs(fprod) : -fabs(fprod)) > 0) gradSign = 1;
  else gradSign = -1;

  // now perform product rule:
  // if F = f_1 * f_2 * .. * f_n, then
  // \nabla F = \sum_i ((\nabla f_i) \prod_{j \neq i} f_j)
  MxDimVector<double, DIM> grad(0), tmp;
  for (int i = 0; i < numShapes; i++) {
    tmp = shapes[i]->gradFunc(p);
    for (int j = 0; j < numShapes; j++)
      if (j != i) tmp *= fVals[j];
    grad += tmp;
  }

  // boy this is an expensive function
  return gradSign * res;
}

#endif
