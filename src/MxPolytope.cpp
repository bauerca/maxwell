
#include "MxPolytope.hpp"

#include <limits>

template<size_t DIM>
MxPolytopeState MxPolytope<DIM>::getState(const MxShape<DIM> & aShape,
const std::vector<MxDimVector<double, DIM> > & verts, std::vector<double> & fVals) const {
  double eps = std::numeric_limits<double>::epsilon();

  fVals.resize(numVerts);

  bool inOn = true, outOn = true, on = true;
  int fs;
  for (size_t i = 0; i < numVerts; ++i) {
    fVals[i] = aShape.func(verts[i]);
    fs = MxUtil::sign(fVals[i]);
    //std::cout << fVals[i] << ", ";
    //if (fVals[i] < -eps) {
    if (fs == -1) {
      inOn = false;
      on = false;
    //} else if (fVals[i] > eps) {
    } else if (fs == 1) {
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
  // polytope, determine if the edge is intersected by the shape by evaluating the shape
  // function at the edge endpoints. If both function values are of different sign,
  // the edge has an intersection point. Find it and store it in edgeXPtsPtrs.
  size_t vertIndx1, vertIndx2, nodeXVertIndx;
  bool hasNodeX;
  double f1, f2;
  int f1s, f2s;
  for (size_t edgeIndx = 0; edgeIndx < numEdges; edgeIndx++) {
    hasNodeX = false;
    vertIndx1 = edgeVertsLists[edgeIndx][0];
    vertIndx2 = edgeVertsLists[edgeIndx][1];
    f1 = fVals[vertIndx1];
    f2 = fVals[vertIndx2];
    f1s = MxUtil::sign(f1);
    f2s = MxUtil::sign(f2);

    if ((f1s == -1 and f2s == 1) or (f1s == 1 and f2s == -1)) {
      edgeX[edgeIndx] = new MxDimVector<double, DIM>(MxUtil::rootFind<MxShape<DIM>, DIM>(aShape, verts[vertIndx1], verts[vertIndx2], 1.e-12, 20));
      edgeXInds.push_back(edgeIndx);
    } else if (f1s == 0 and f2s == 1) {
      //hasNodeX = true;
      //nodeXVertIndx = vertIndx1;
      edgeX[edgeIndx] = new MxDimVector<double, DIM>(verts[vertIndx1]);
      edgeXInds.push_back(edgeIndx);
    } else if (f2s == 0 and f1s == 1) {
      //hasNodeX = true;
      //nodeXVertIndx = vertIndx2;
      edgeX[edgeIndx] = new MxDimVector<double, DIM>(verts[vertIndx2]);
      edgeXInds.push_back(edgeIndx);
    }

/*
    if (hasNodeX) {
      typename std::vector<MxDimVector<double, DIM> *>::const_iterator iter;
      bool insert = true;
      for (size_t i = 0; i < edgeXInds.size(); ++i) {
        if (*edgeX[edgeXInds[i]] == verts[nodeXVertIndx]) {
          insert = false;
        }
      }
      if (insert) {
        edgeX[edgeIndx] = new MxDimVector<double, DIM>(verts[nodeXVertIndx]);
        edgeXInds.push_back(edgeIndx);
      }
    }
*/
  }
}


template class MxPolytope<1>;
template class MxPolytope<2>;
template class MxPolytope<3>;
