
#include <cmath>
#include <utility>
#include <limits>

#include "MxConvexPolyhedron.h"
#include "MxUtil.hpp"

#include "MxShape.hpp"

//MxConvexPolyhedron::MxConvexPolyhedron(std::vector<MxDimVector<double, 3> > vertices,
//MxDynMatrix<size_t> vertEdgeConn,
//MxDynMatrix<size_t> 
//double sideLength1, MxDimVector<double, 3> direction1, 



//MxConvexPolyhedron::MxConvexPolyhedron(double sideLength1, MxDimVector<double, 3> direction1, 
//double sideLength2, MxDimVector<double, 3> direction2,
//double sideLength3, MxDimVector<double, 3> direction3, MxDimVector<double, 3> midPoint) : 
//numVerts(4), l1(sideLength1), l2(sideLength2), l3(sideLength3), 
//d1(direction1), d2(direction2), d3(direction3), midPt(midPoint) {
//
//  vertices = this->getVertices(midPoint);
//   
//  edges.resize(12);
//  // length1 edges
//  edges[0] = std::make_pair(0, 4);
//  edges[1] = std::make_pair(1, 5);
//  edges[2] = std::make_pair(2, 6);
//  edges[3] = std::make_pair(3, 7);
//  // length2 edges
//  edges[4] = std::make_pair(0, 2);
//  edges[5] = std::make_pair(1, 3);
//  edges[6] = std::make_pair(4, 6);
//  edges[7] = std::make_pair(5, 7);
//  // length3 edges
//  edges[8] = std::make_pair(0, 1);
//  edges[9] = std::make_pair(2, 3);
//  edges[10] = std::make_pair(4, 5);
//  edges[11] = std::make_pair(6, 7);
//
//  faces.resize(6);
//  // from length1, length2 edges
//  faces[0] = 
//}

void MxConvexPolyhedron::setLists() {
  this->setEdgeFacesLists();
  this->setFaceEdgesLists();
  this->setEdgeVertsLists();
  this->setFaceVertsLists();
}

//std::vector<size_t> MxConvexPolyhedron::getFaceEdgeIndxs(size_t faceIndx) const {
//  std::vector<size_t> res;
//  for (size_t edgeIndx = 0; edgeIndx < numEdges; edgeIndx++)
//    if (edgeFaceMatrix(edgeIndx, faceIndx) != 0)
//      res.push_back(edgeIndx);
//  return res;
//}
//
//std::vector<size_t> MxConvexPolyhedron::getEdgeFaceIndxs(size_t edgeIndx) const {
//  std::vector<size_t> res;
//  for (size_t faceIndx = 0; faceIndx < numFaces; faceIndx++)
//    if (edgeFaceMatrix(edgeIndx, faceIndx) != 0)
//      res.push_back(faceIndx);
//  return res;
//}
//
//std::vector<size_t> MxConvexPolyhedron::getFaceVertIndxs(size_t faceIndx) const {
//  std::vector<size_t> res;
//  for (size_t vertIndx = 0; vertIndx < numVerts; vertIndx++)
//    if (faceVertMatrix(faceIndx, vertIndx) != 0)
//      res.push_back(vertIndx);
//  return res;
//}
//
//std::vector<size_t> MxConvexPolyhedron::getEdgeVertIndxs(size_t edgeIndx) const {
//  std::vector<size_t> res;
//  for (size_t vertIndx = 0; vertIndx < numVerts; vertIndx++)
//    if (faceVertMatrix(faceIndx, vertIndx) != 0)
//      res.push_back(vertIndx);
//  return res;
//}

void MxConvexPolyhedron::setEdgeVertsLists() {
  edgeVertsLists.resize(numEdges);
  for (size_t edgeIndx = 0; edgeIndx < numEdges; edgeIndx++) {
    std::vector<size_t> edgeVerts;
    for (size_t vertIndx = 0; vertIndx < this->numVerts; vertIndx++)
      if (edgeVertMatrix(edgeIndx, vertIndx) != 0)
        edgeVerts.push_back(vertIndx);
    edgeVertsLists[edgeIndx] = edgeVerts;
  }
}

void MxConvexPolyhedron::setEdgeFacesLists() {
  edgeFacesLists.resize(numEdges);
  for (size_t edgeIndx = 0; edgeIndx < numEdges; edgeIndx++) {
    std::vector<size_t> edgeFaces;
    for (size_t faceIndx = 0; faceIndx < numFaces; faceIndx++)
      if (edgeFaceMatrix(edgeIndx, faceIndx) != 0)
        edgeFaces.push_back(faceIndx);
    edgeFacesLists[edgeIndx] = edgeFaces;
  }
}

void MxConvexPolyhedron::setFaceEdgesLists() {
  faceEdgesLists.resize(numFaces);
  for (size_t faceIndx = 0; faceIndx < numFaces; faceIndx++) {
    std::vector<size_t> faceEdges;
    for (size_t edgeIndx = 0; edgeIndx < numEdges; edgeIndx++)
      if (edgeFaceMatrix(edgeIndx, faceIndx) != 0)
        faceEdges.push_back(edgeIndx);
    faceEdgesLists[faceIndx] = faceEdges;
  }
}

void MxConvexPolyhedron::setFaceVertsLists() {
  faceVertsLists.resize(numFaces);
  for (size_t faceIndx = 0; faceIndx < numFaces; faceIndx++) {
    std::vector<size_t> faceVerts;
    for (size_t vertIndx = 0; vertIndx < this->numVerts; vertIndx++)
      if (faceVertMatrix(faceIndx, vertIndx) != 0)
        faceVerts.push_back(vertIndx);
    faceVertsLists[faceIndx] = faceVerts;
  }
}

double MxConvexPolyhedron::volumeFraction(const MxShape<3> & aShape, const MxDimVector<double, 3> & p) const {

  double res;

  std::vector<MxDimVector<double, 3> > verts = this->getVertices(p);

  std::vector<double> fVals;
  // fill fVals
  MxPolytopeState myState = this->getState(aShape, verts, fVals);
  //std::cout << "fVals: "; MxUtil::printStdVector(fVals);

  std::vector<MxDimVector<double, 3> *> edgeX;
  std::vector<size_t> edgeXInds;

  switch (myState) {
    case ON:
      if (aShape.func(p) > MxUtil::dEps) res = 1;
      else res = 0;
      break;
    case IN:
      res = 1;
      break;
    case OUT:
      res = 0;
      break;
    case CUT1:
      this->getEdgeX(aShape, verts, fVals, edgeXInds, edgeX);
      res = this->volumeFraction(aShape, verts, fVals, edgeXInds, edgeX);
      break;
  }

  //std::cout << "volume fraction for polyhedron: " << res << "\n";

  std::vector<MxDimVector<double, 3> *>::iterator iter;
  for (iter = edgeX.begin(); iter != edgeX.end(); ++iter)
    delete *iter;

  return res;
}

double MxConvexPolyhedron::volumeFraction(
const MxShape<3> & aShape,
const std::vector<MxDimVector<double, 3> > & verts, 
const std::vector<double> & fVals, 
const std::vector<size_t> & edgeXInds,
const std::vector<MxDimVector<double, 3> *> & edgeX) const {

  // This will be the final volume of the portion of the convex polyhedron that is inside
  // the shape. The volume will be calculated as a sum of pyramidal volumes.
  double vol = 0;

  MxDimVector<double, 3> w1, w2; // work vectors

  // the following block is for the calculation of face area vectors. For every face,
  // we first determine if any part of the face is inside the given shape. If so, we 
  // then determine whether that face is intersected by the boundary of the shape. If
  // the face is intersected, 

  std::vector<MxDimVector<double, 3> *> faceAreaVecPtrs(numFaces, 0);
  std::vector<const MxDimVector<double, 3> *> faceVertPtrs(numFaces, 0);
  std::vector<double> faceCutLengths(numFaces, 0);
  const MxDimVector<double, 3> * v0 = 0;
  const MxDimVector<double, 3> * v1 = 0;
  const MxDimVector<double, 3> * v2 = 0;
  MxDimVector<double, 3> * areaVec = 0;

  bool inside, cut;
  std::vector<size_t> faceEdgeIndxs, faceVertIndxs, edgeFaceIndxs;
  std::vector<size_t>::iterator vertIndxIter, edgeIndxIter;

  size_t vertIndx1, vertIndx2;

  // loop over all faces
  for (size_t faceIndx = 0; faceIndx < numFaces; faceIndx++) {
    faceEdgeIndxs = faceEdgesLists[faceIndx];
    faceVertIndxs = faceVertsLists[faceIndx];

    //determine if any part of face is inside
    inside = false;
    for (vertIndxIter = faceVertIndxs.begin(); vertIndxIter != faceVertIndxs.end(); vertIndxIter++)
      if (fVals[*vertIndxIter] > 0) {inside = true; break;}

    // go to next face early if current face is completely outside shape
    if (not inside) continue;

    // determine if face is cut. If so, get a handle for the edge
    // intersection point to be used as the base vertex for all triangular
    // areas to be calculated for this face.
    //
    cut = false;
    for (edgeIndxIter = faceEdgeIndxs.begin(); edgeIndxIter != faceEdgeIndxs.end(); edgeIndxIter++) {
      if (edgeX[*edgeIndxIter] != 0) {
        cut = true; 
        v0 = edgeX[*edgeIndxIter];
        break;
      }
    }

    //std::cout << "Calculating area vector for face " << faceIndx << "\n";

    areaVec = new MxDimVector<double, 3>(0);
    faceAreaVecPtrs[faceIndx] = areaVec;

    if (cut) {
      // grab a sample vertex from this face
      faceVertPtrs[faceIndx] = v0;
      // current face is cut. Loop over all face edges, using only those
      // that are cut or fully inside
      for (edgeIndxIter = faceEdgeIndxs.begin(); edgeIndxIter != faceEdgeIndxs.end(); edgeIndxIter++) {
        vertIndx1 = edgeVertsLists[*edgeIndxIter][0];
        vertIndx2 = edgeVertsLists[*edgeIndxIter][1];
        w1 = verts[vertIndx1];
        w2 = verts[vertIndx2];

        //std::cout << edgeX[*edgeIndxIter] << std::endl;
        // if edge is cut and different
        if ((edgeX[*edgeIndxIter] != 0) && (edgeX[*edgeIndxIter] != v0)) {
          v1 = edgeX[*edgeIndxIter];
          v2 = fVals[vertIndx1] > 0 ? &w1 : &w2;

          faceCutLengths[faceIndx] = (*v0 - *v1).norm();
        }
        // if edge is not cut and inside
        else if ((fVals[vertIndx1] > MxUtil::dEps) or (fVals[vertIndx2] > MxUtil::dEps)) {
          v1 = &w1;
          v2 = &w2;
        } 
        // otherwise edge must be outside, so don't calculate a triangle area
        else continue;
      
        //std::cout << "triangle vertices:\n";
        //v0->print();
        //v1->print();
        //v2->print();
        // If we get here, is current triangle area vec pointing the same direction?
        MxDimVector<double, 3> triAreaVec(cross(*v1 - *v0, *v2 - *v0));
        *areaVec += areaVec->dot(triAreaVec) > 0 ? triAreaVec : -triAreaVec;
      }
    }

    // if face is fully inside boundary
    else {
      // any 'ol face vertex will do
      size_t baseVertIndx = faceVertIndxs[0];
      v0 = &verts[baseVertIndx];
      faceVertPtrs[faceIndx] = v0;
      // loop over all face edges
      for (edgeIndxIter = faceEdgeIndxs.begin(); edgeIndxIter != faceEdgeIndxs.end(); edgeIndxIter++) {
        // get edge endpoints
        vertIndx1 = edgeVertsLists[*edgeIndxIter][0];
        vertIndx2 = edgeVertsLists[*edgeIndxIter][1];
        w1 = verts[vertIndx1];
        w2 = verts[vertIndx2];

        if ((vertIndx1 != baseVertIndx) and (vertIndx2 != baseVertIndx)) {
          v1 = &w1;
          v2 = &w2;
        }
        else continue;

        //std::cout << "triangle vertices:\n";
        //v0->print();
        //v1->print();
        //v2->print();

        //is current triangle area vec pointing the same direction?
        MxDimVector<double, 3> triAreaVec(cross(*v1 - *v0, *v2 - *v0));
        *areaVec += areaVec->dot(triAreaVec) > 0 ? triAreaVec : -triAreaVec;
      }
    }

    *areaVec *= 0.5;
  }

  // have intersection points, weights. Need area inside shape for each face
  // For each edge intersection point, we calculate the volume of the pyramid 
  // formed by connecting the portion of each face inside the shape with
  // the edge intersection point.
  //
  // sum up all tetrahedra
  double wt, wtSum = 0;
  for (size_t edgeIndx = 0; edgeIndx < numEdges; edgeIndx++) {
    if (edgeX[edgeIndx] == 0) continue;

    edgeFaceIndxs = edgeFacesLists[edgeIndx];
    wt = faceCutLengths[edgeFaceIndxs[0]] + faceCutLengths[edgeFaceIndxs[1]];
    wtSum += wt;

    for (size_t faceIndx = 0; faceIndx < numFaces; faceIndx++) {
      if (faceAreaVecPtrs[faceIndx] != 0) {
        //faceAreaVecPtrs[faceIndx]->print();
        vol += wt * fabs((*faceVertPtrs[faceIndx] - *edgeX[edgeIndx]).dot(*faceAreaVecPtrs[faceIndx])) / 3.;
      }
    }
  }

  // delete stuff
  for (size_t faceIndx = 0; faceIndx < numFaces; faceIndx++)
    delete faceAreaVecPtrs[faceIndx];

  vol /= wtSum;

  return vol / this->volume();
}


