
#include "MxSegment.h"
#include "MxUtil.hpp"
#include "MxShape.hpp"
#include "MxDynMatrix.hpp"


template<size_t DIM>
MxSegment<DIM>::MxSegment(const MxDimVector<double, DIM> & p1, const MxDimVector<double, DIM> & p2) : 
numVerts(2) {
  len = (p2 - p1).norm();
  dir = (p2 - p1) / len;
  edgeVertMatrix = MxDynMatrix<size_t>(1, 2, 1); // = [1, 1]
}

template<size_t DIM>
MxSegment<DIM>::MxSegment(double aLen, const MxDimVector<double, DIM> & aDir) : 
numVerts(2), len(aLen), dir(aDir / aDir.norm()), edgeVertMatrix(MxDynMatrix<size_t>(1, 2, 1)) {}


template<size_t DIM>
double MxSegment<DIM>::volumeFraction(const MxShape<DIM> & aShape, const MxDimVector<double, DIM> & midPt) const {

  MxDimVector<double, DIM> p1(midPt + 0.5 * len * dir);
  MxDimVector<double, DIM> p2(midPt - 0.5 * len * dir);
  double f1 = aShape.func(p1);
  double f2 = aShape.func(p2);

  if      ((f1 > MxUtil::dEps and f2 > -MxUtil::dEps) or 
           (f2 > MxUtil::dEps and f1 > -MxUtil::dEps)) return 1.0;     // one in, other in or on
  else if ((f1 < -MxUtil::dEps and f2 < MxUtil::dEps) or
           (f2 < -MxUtil::dEps and f1 < MxUtil::dEps)) return 0.0;     // one out, other out or on
  else if (fabs(f1) < MxUtil::dEps and fabs(f2) < MxUtil::dEps) {      // both on
    if (aShape.func(midPt) > MxUtil::dEps) return 1.0;
    else return 0.0;
  }
  else {                                                                         // bracketing
    MxDimVector<double, DIM> p(MxUtil::rootFind<MxShape<DIM>, DIM>(aShape, p1, p2, 1.e-12, 20));
    if (f1 > MxUtil::dEps)
      return (p - p1).norm() / len;
    else
      return (p - p2).norm() / len;
  }

#if 0
  if (f1 * f2 > 0.0) {
    if (f1 > 0.0)
      return 1.0;
    else
      return 0.0;
  }
  else if (f1 == 0.0 and f2 == 0.0) {
    if (aShape.func(midPt) > 0.0) return 1.0;
    else return 0.0;
  }
  else if (f1 * f2 == 0.0) {
    if (f1 > 0.0 or f2 > 0.0) return 1.0;
    else return 0.0;
  }
  else {
    MxDimVector<double, DIM> p(MxUtil::rootFind<MxShape<DIM>, DIM>(aShape, p1, p2, 1.e-12, 20));
    if (f1 > 0.0)
      return (p - p1).norm() / len;
    else
      return (p - p2).norm() / len;
  }
#endif

}

// this function assumes that there is an edge intersection
template<size_t DIM>
double MxSegment<DIM>::volumeFraction(const MxShape<DIM> & aShape, 
const std::vector<MxDimVector<double, DIM> > & verts,
const std::vector<double> & fVals,
const std::vector<size_t> & edgeXInds,
const std::vector<MxDimVector<double, DIM> *> & edgeX) const {

  if (fVals[0] > MxUtil::dEps)
    return (*edgeX[0] - verts[0]).norm() / len;
  else
    return (*edgeX[0] - verts[1]).norm() / len;

}

template<size_t DIM>
std::vector<MxDimVector<double, DIM> > MxSegment<DIM>::getVertices(const MxDimVector<double, DIM> & midPt) const {
  std::vector<MxDimVector<double, DIM> > res;
  res.push_back(midPt - 0.5 * len * dir);
  res.push_back(midPt + 0.5 * len * dir);
  return res;
}

template<size_t DIM>
void MxSegment<DIM>::getEdgeX(const MxShape<DIM> & aShape,
const std::vector<MxDimVector<double, DIM> > & verts,
const std::vector<double> & fVals,
std::vector<size_t> & edgeXInds,
std::vector<MxDimVector<double, DIM> *> & edgeX) const {
  edgeXInds.resize(1);
  edgeX.resize(1, 0);
  edgeXInds.push_back(0);
  edgeX.push_back(new MxDimVector<double, DIM>(MxUtil::rootFind<MxShape<DIM>, DIM>(aShape, verts[0], verts[1], 1.e-12, 20)));
}

// instantiation

template class MxSegment<1>;
template class MxSegment<2>;
template class MxSegment<3>;
