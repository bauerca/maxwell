
#include "MxCartBox.hpp"


MxCartBox::MxCartBox(double aLx, double aLy, double aLz) : 
lx(aLx), ly(aLy), lz(aLz) {
  numVerts = 8;
  numEdges = 12;
  numFaces = 6;
  setConnMatrices();
}

MxCartBox::MxCartBox(const MxDimVector<double, 3> & l) : 
lx(l[0]), ly(l[1]), lz(l[2]) {
  numVerts = 8;
  numEdges = 12;
  numFaces = 6;
  setConnMatrices();
}

void MxCartBox::setConnMatrices() {
  edgeVertMatrix = MxDynMatrix<size_t>(numEdges, numVerts, 0);
  // x-edges
  edgeVertMatrix(0, 0) = 1; edgeVertMatrix(0, 4) = 1;
  edgeVertMatrix(1, 1) = 1; edgeVertMatrix(1, 5) = 1;
  edgeVertMatrix(2, 2) = 1; edgeVertMatrix(2, 6) = 1;
  edgeVertMatrix(3, 3) = 1; edgeVertMatrix(3, 7) = 1;
  // y-edges
  edgeVertMatrix(4, 0) = 1; edgeVertMatrix(4, 2) = 1;
  edgeVertMatrix(5, 1) = 1; edgeVertMatrix(5, 3) = 1;
  edgeVertMatrix(6, 4) = 1; edgeVertMatrix(6, 6) = 1;
  edgeVertMatrix(7, 5) = 1; edgeVertMatrix(7, 7) = 1;
  // z-edges
  edgeVertMatrix(8, 0) = 1; edgeVertMatrix(8, 1) = 1;
  edgeVertMatrix(9, 2) = 1; edgeVertMatrix(9, 3) = 1;
  edgeVertMatrix(10, 4) = 1; edgeVertMatrix(10, 5) = 1;
  edgeVertMatrix(11, 6) = 1; edgeVertMatrix(11, 7) = 1;

  edgeFaceMatrix = MxDynMatrix<size_t>(numEdges, numFaces, 0);
  // x-edges                                          // edge  : face
  edgeFaceMatrix(0, 0) = 1; edgeFaceMatrix(0, 2) = 1; // -y,-z : -z ; -y
  edgeFaceMatrix(1, 1) = 1; edgeFaceMatrix(1, 2) = 1; // -y,+z : +z ; -y
  edgeFaceMatrix(2, 0) = 1; edgeFaceMatrix(2, 3) = 1; // +y,-z : -z ; +y
  edgeFaceMatrix(3, 1) = 1; edgeFaceMatrix(3, 3) = 1; // +y,+z : +z ; +y
  // y-edges
  edgeFaceMatrix(4, 0) = 1; edgeFaceMatrix(4, 4) = 1; // -x,-z : -z ; -x
  edgeFaceMatrix(5, 1) = 1; edgeFaceMatrix(5, 4) = 1; // -x,+z : +z ; -x
  edgeFaceMatrix(6, 0) = 1; edgeFaceMatrix(6, 5) = 1; // +x,-z : -z ; +x
  edgeFaceMatrix(7, 1) = 1; edgeFaceMatrix(7, 5) = 1; // +x,+z : +z ; +x
  // z-edges
  edgeFaceMatrix(8, 2) = 1; edgeFaceMatrix(8, 4) = 1;   // -x,-y : -y ; -x
  edgeFaceMatrix(9, 3) = 1; edgeFaceMatrix(9, 4) = 1;   // -x,+y : +y ; -x
  edgeFaceMatrix(10, 2) = 1; edgeFaceMatrix(10, 5) = 1; // +x,-y : -y ; +x
  edgeFaceMatrix(11, 3) = 1; edgeFaceMatrix(11, 5) = 1; // +x,+y : +y ; +x

  //edgeFaceMatrix.print();

  faceVertMatrix = (edgeFaceMatrix.transpose()) * (edgeVertMatrix);

  this->setLists();

}

