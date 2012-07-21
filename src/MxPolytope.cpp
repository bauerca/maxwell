
#include "MxPolytope.hpp"

#include <limits>

template<size_t DIM>
MxPolytopeState MxPolytope<DIM>::getState(const MxShape<DIM> & aShape,
const std::vector<MxDimVector<double, DIM> > & verts, std::vector<double> & fVals) const {
  double eps = std::numeric_limits<double>::epsilon();

  fVals.resize(numVerts);

  bool inOn = true, outOn = true, on = true;
  for (size_t i = 0; i < numVerts; ++i) {
    fVals[i] = aShape.func(verts[i]);
    //std::cout << fVals[i] << ", ";
    if (fVals[i] < -eps) {
      inOn = false;
      on = false;
    }
    if (fVals[i] > eps) {
      outOn = false;
      on = false;
    }
  }
  //std::cout << "\n";

  if (on) return ON;
  else if (inOn) return IN;
  else if (outOn) return OUT;
  else return CUT1;
}



template<size_t DIM>
void MxPolytope<DIM>::getEdgeX(
const MxShape<DIM> & aShape, 
const std::vector<MxDimVector<double, DIM> > & verts, 
const std::vector<double> & fVals, 
std::vector<size_t> & edgeXInds,
std::vector<MxDimVector<double, DIM> *> & edgeX) const {
  edgeXInds.clear();
  edgeX.resize(numEdges, 0);

  // First, we enumerate the edge intersection points. For every edge of the 
  // parallelogram, determine if the edge is intersected by the shape by evaluating the shape
  // function at the edge endpoints. If both function values are of different sign,
  // the edge has an intersection point. Find it and store it in edgeXPtsPtrs.
  size_t vertIndx1, vertIndx2;
  double f1, f2;
  for (size_t edgeIndx = 0; edgeIndx < numEdges; edgeIndx++) {
    vertIndx1 = edgeVertsLists[edgeIndx][0];
    vertIndx2 = edgeVertsLists[edgeIndx][1];
    f1 = fVals[vertIndx1];
    f2 = fVals[vertIndx2];
    if (((f1 < MxUtil::dEps) and (f2 > MxUtil::dEps)) or ((f1 > MxUtil::dEps) and (f2 < MxUtil::dEps))) {
      if (abs(f1) < MxUtil::dEps)
        edgeX[edgeIndx] = new MxDimVector<double, DIM>(verts[vertIndx1]);
      else if (abs(f2) < MxUtil::dEps)
        edgeX[edgeIndx] = new MxDimVector<double, DIM>(verts[vertIndx2]);
      else
        edgeX[edgeIndx] = new MxDimVector<double, DIM>(MxUtil::rootFind<MxShape<DIM>, DIM>(aShape, verts[vertIndx1], verts[vertIndx2], 1.e-12, 20));
      edgeXInds.push_back(edgeIndx);
    }
  }
}


template class MxPolytope<1>;
template class MxPolytope<2>;
template class MxPolytope<3>;
