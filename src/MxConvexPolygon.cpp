
#include <cmath>

#include "MxConvexPolygon.h"
#include "MxUtil.hpp"

#include "MxShape.hpp"
//#include "MxShapeDefs.h"
#include "MxHalfSpace.hpp"


template<size_t DIM>
void MxConvexPolygon<DIM>::setEdgeVertsLists() {
  this->edgeVertsLists.resize(this->numEdges);
  for (size_t edgeIndx = 0; edgeIndx < this->numEdges; edgeIndx++) {
    std::vector<size_t> edgeVerts;
    for (size_t vertIndx = 0; vertIndx < this->numVerts; vertIndx++)
      if (edgeVertMatrix(edgeIndx, vertIndx) != 0)
        edgeVerts.push_back(vertIndx);
    this->edgeVertsLists[edgeIndx] = edgeVerts;
  }
}

template<size_t DIM>
double MxConvexPolygon<DIM>::volumeFraction(
const MxShape<DIM> & aShape,
const std::vector<MxDimVector<double, DIM> > & verts, 
const std::vector<double> & fVals, 
const std::vector<size_t> & edgeXInds,
const std::vector<MxDimVector<double, DIM> *> & edgeX) const {
  double area = 0.0;

  MxDimVector<double, DIM> w1, w2;

  // Each triangular area needs 3 vertices. 
  MxDimVector<double, DIM> v0(this->getCorner(aShape, verts, edgeXInds, edgeX));
  MxDimVector<double, DIM> * v1, * v2;

  // Now loop over all edges again. For each edge, we calculate a trianglar area.
  // If the edge is not intersected, then the endpoints of the edge are the 
  // remaining two vertices for the triangle. If the edge is intersected, then
  // the endpoint that is outside the shape is thrown away
  size_t vertIndx1, vertIndx2;
  int f1s, f2s;
  for (size_t edgeIndx = 0; edgeIndx < this->numEdges; edgeIndx++) {
    vertIndx1 = this->edgeVertsLists[edgeIndx][0];
    vertIndx2 = this->edgeVertsLists[edgeIndx][1];
    w1 = verts[vertIndx1];
    w2 = verts[vertIndx2];

    f1s = MxUtil::sign(fVals[vertIndx1]);
    f2s = MxUtil::sign(fVals[vertIndx2]);

    // if edge is cut,
    // use endpoint inside shape and intersection point
    if (edgeX[edgeIndx] != 0) {
      v1 = f1s != -1 ? &w1 : edgeX[edgeIndx];
      v2 = f2s != -1 ? &w2 : edgeX[edgeIndx];
      //v1 = fVals[vertIndx1] > 0 ? &w1 : edgeX[edgeIndx];
      //v2 = fVals[vertIndx2] > 0 ? &w2 : edgeX[edgeIndx];
    }
    // if edge is inside use both edge endpoints
    //else if ((fVals[vertIndx1] > MxUtil::dEps) or (fVals[vertIndx2] > MxUtil::dEps)) {
    else if (f1s == 1 or f2s == 1) {
      v1 = &w1;
      v2 = &w2;
    }
    // if edge is outside (or on?), skip it!
    else continue;

    // Add the area for this triangle
    //   area of triangle = 0.5 * |(v1 - v0) x (v2 - v0)|
    area += 0.5 * cross(*v1 - v0, *v2 - v0).norm();
  }

  return area / this->volume();
}

template<>
MxDimVector<double, 2> MxConvexPolygon<2>::getCorner(
const MxShape<2> & aShape, 
const std::vector<MxDimVector<double, 2> > & verts,
const std::vector<size_t> & edgeXInds,
const std::vector<MxDimVector<double, 2> *> & edgeX) const {
  MxDimVector<double, 2> tmp;

  std::vector<MxDimVector<double, 2> > edgeXNormals;
  tmp = aShape.gradFunc(*edgeX[edgeXInds[0]]);
  tmp /= tmp.norm();
  edgeXNormals.push_back(tmp);
  tmp = aShape.gradFunc(*edgeX[edgeXInds[1]]);
  tmp /= tmp.norm();
  edgeXNormals.push_back(tmp);

  // first check angle between normals at intersection points
  double cosTheta = edgeXNormals[0].dot(edgeXNormals[1]);
    
  // if angle is too large (should user define this angle?), find the 'corner' (could just
  // be a really curvy surface)
  //if (cosTheta < 1.0 / sqrt(2.0)) {
  if (cosTheta < cos(0.9 * 0.5 * MxUtil::pi)) {
    std::cout << "n1 = "; edgeXNormals[0].print(); std::cout << "  at point: "; edgeX[edgeXInds[0]]->print();
    std::cout << "n2 = "; edgeXNormals[1].print(); std::cout << "  at point: "; edgeX[edgeXInds[1]]->print();
    std::cout << "n1 . n2 = " << cosTheta << "\n";
    MxDimMatrix<double, 2> mat, matInv;
    MxDimVector<double, 2> b;

    mat.setRow(0, edgeXNormals[0]);
    mat.setRow(1, edgeXNormals[1]);
    //std::cout << "matrix to invert to find corner\n";
    //mat.print();
    matInv = mat.inv();
    
    b[0] = edgeXNormals[0].dot(*edgeX[edgeXInds[0]]);
    b[1] = edgeXNormals[1].dot(*edgeX[edgeXInds[1]]);

    tmp = matInv * b;
    std::cout << "corner guess:"; tmp.print();

    // check that guess is within bounds of polygon

// 0 to use corner guess
#if 0
    //std::cout << "vertices:\n";
    //for (size_t i = 0; i < verts.size(); ++i) {
    //  std::cout << "  "; verts[i].print();
    //}
    std::cout << "corner guess:"; tmp.print();
    std::cout << "corner guess func val: " << aShape.func(tmp) << "\n";

    // use 'corner' as initial guess to multidimensional root find
    int shpIndx1, shpIndx2;
    aShape.func(*edgeX[edgeXInds[0]], shpIndx1);
    aShape.func(*edgeX[edgeXInds[1]], shpIndx2);
    MxShape<2> * shp1 = aShape.getSubShape(shpIndx1);
    MxShape<2> * shp2 = aShape.getSubShape(shpIndx2);
    std::cout << "Shape1 name: " << shp1->getName() << "; index: " << shpIndx1 << "\n";
    std::cout << "Shape2 name: " << shp2->getName() << "; index: " << shpIndx2 << "\n";

    std::vector<const MxShape<2> *> shapes;
    shapes.push_back(shp1);
    shapes.push_back(shp2);
    tmp = MxUtil::ndRootFind<MxShape<2>, 2>(shapes, tmp, 1.e-12, 20);

    std::cout << "corner exact:"; tmp.print();
    if (fabs(aShape.func(tmp)) > 1.e-10) {
      std::cout << "corner exact func val: " << aShape.func(tmp) << "\n";
      throw 1;
    }
      
    std::cout << std::endl;

    delete shp1;
    delete shp2;
#endif

    return tmp;
  }
  else {
    return *edgeX[edgeXInds[0]];
  }

}

template<>
MxDimVector<double, 3> MxConvexPolygon<3>::getCorner(
const MxShape<3> & aShape, 
const std::vector<MxDimVector<double, 3> > & verts,
const std::vector<size_t> & edgeXInds,
const std::vector<MxDimVector<double, 3> *> & edgeX) const {
  MxDimVector<double, 3> tmp;

  //std::cout << "edgeXInds size=" << edgeXInds.size() << "\n";
  //for (int i = 0; i < edgeXInds.size(); ++i) {
  //  std::cout << "  edgeXInd " << i << ": " << edgeXInds[i] << "\n";
  //  std::cout << "  edgeX " << i << ": " << *edgeX[edgeXInds[i]] << "\n";
  //}

  std::vector<MxDimVector<double, 3> > edgeXNormals;
  tmp = aShape.gradFunc(*edgeX[edgeXInds[0]]);
  tmp /= tmp.norm();
  edgeXNormals.push_back(tmp);
  tmp = aShape.gradFunc(*edgeX[edgeXInds[1]]);
  tmp /= tmp.norm();
  edgeXNormals.push_back(tmp);

  // first check angle between normals at intersection points
  double cosTheta = edgeXNormals[0].dot(edgeXNormals[1]);
    
  // if angle is too large (should user define this angle?), find the 'corner' (could just
  // be a really curvy surface)
  //if (cosTheta < 1.0 / sqrt(2.0)) {
  //if (cosTheta < cos(0.9 * 0.5 * MxUtil::pi)) {
  if (false) {
    MxVecD3 polyNorm(cross(verts[1] - verts[0], verts[2] - verts[0]));

    std::cout << "n1 = "; edgeXNormals[0].print(); std::cout << "  at point: "; edgeX[edgeXInds[0]]->print();
    std::cout << "n2 = "; edgeXNormals[1].print(); std::cout << "  at point: "; edgeX[edgeXInds[1]]->print();
    std::cout << "n1 . n2 = " << cosTheta << "\n";
    MxDimMatrix<double, 3> mat, matInv;
    MxDimVector<double, 3> b;

    mat.setRow(0, edgeXNormals[0]);
    mat.setRow(1, edgeXNormals[1]);
    mat.setRow(2, polyNorm);
    //std::cout << "matrix to invert to find corner\n";
    //mat.print();
    matInv = mat.inv();
    
    b[0] = edgeXNormals[0].dot(*edgeX[edgeXInds[0]]);
    b[1] = edgeXNormals[1].dot(*edgeX[edgeXInds[1]]);
    b[2] = polyNorm.dot(*edgeX[edgeXInds[0]]); // any point on the polygon will do

    tmp = matInv * b;
    //std::cout << "vertices:\n";
    //for (size_t i = 0; i < verts.size(); ++i) {
    //  std::cout << "  "; verts[i].print();
    //}
    std::cout << "corner guess:"; tmp.print();
    std::cout << "corner guess func val: " << aShape.func(tmp) << "\n";

// set 0 to use corner guess
#if 0
    // use 'corner' as initial guess to multidimensional root find
    int shpIndx1, shpIndx2;
    aShape.func(*edgeX[edgeXInds[0]], shpIndx1);
    aShape.func(*edgeX[edgeXInds[1]], shpIndx2);
    MxHalfSpace<3> polyPlane(verts[0], polyNorm);
    MxShape<3> * shp1 = aShape.getSubShape(shpIndx1);
    MxShape<3> * shp2 = aShape.getSubShape(shpIndx2);
    std::cout << "Shape1 name: " << shp1->getName() << "; index: " << shpIndx1 << "\n";
    std::cout << "Shape2 name: " << shp2->getName() << "; index: " << shpIndx2 << "\n";

    std::vector<const MxShape<3> *> shapes;
    shapes.push_back(&polyPlane);
    shapes.push_back(shp1);
    shapes.push_back(shp2);
    tmp = MxUtil::ndRootFind<MxShape<3>, 3>(shapes, tmp, 1.e-12, 20);

    std::cout << "corner exact:"; tmp.print();
    std::cout << std::endl;

    delete shp1;
    delete shp2;
#endif

    return tmp;
  }
  else {
    return *edgeX[edgeXInds[0]];
  }

}


template<size_t DIM>
double MxConvexPolygon<DIM>::volumeFraction(const MxShape<DIM> & aShape, const MxDimVector<double, DIM> & p) const {

  double res;

  std::vector<MxDimVector<double, DIM> > verts = this->getVertices(p);

  std::vector<double> fVals;
  // fill fVals
  MxPolytopeState myState = this->getState(aShape, verts, fVals);

  std::vector<MxDimVector<double, DIM> *> edgeX;
  std::vector<size_t> edgeXInds;

  switch (myState) {
    case ON:
      if (MxUtil::sign(aShape.func(p)) == 1) res = 1;
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

  typename std::vector<MxDimVector<double, DIM> *>::iterator iter;
  for (iter = edgeX.begin(); iter != edgeX.end(); ++iter)
    delete *iter;

  return res;
}


// explicit template instantiation
template class MxConvexPolygon<2>;
template class MxConvexPolygon<3>;

