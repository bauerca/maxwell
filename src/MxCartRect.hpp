#ifndef MX_CART_RECT
#define MX_CART_RECT

#include "MxDimVector.hpp"
#include "MxDynMatrix.hpp"
#include "MxConvexPolygon.h"


template<size_t DIM>
class MxCartRect : public MxConvexPolygon<DIM> {
  public:
    MxCartRect(char aType, double aL1, double aL2);

    // don't include if you want to use the inherited method
    //virtual double volumeFraction(const MxShape<DIM> & aShape, const MxDimVector<double, DIM> & midPt) const;

    virtual double volume() const {return l1 * l2;}

    //virtual MxPolytope<DIM> intersection(const MxShape<DIM> & aShape) const;
    
    virtual std::vector<MxDimVector<double, DIM> > getVertices(const MxDimVector<double, DIM> & p) const;

  private:

    char type;

    double l1, l2;

    MxDimVector<double, DIM> d1, d2;

};

template<size_t DIM>
MxCartRect<DIM>::MxCartRect(char aType, double aL1, double aL2) : 
type(aType), l1(aL1), l2(aL2), d1(0), d2(0) {
  this->numVerts = 4;
  this->numEdges = 4;

  if (DIM == 1) {
    std::cout << "MxCartRect: in 1D, use MxSegment instead.";
    throw 1;
  }
  else if (DIM == 2) {
    if (type != 'z') {
      std::cout << "MxCartRect: Can only have 'z' faces in 2D simulation.";
      throw 1;
    }
    d1[0] = 1; d1[1] = 0;
    d2[0] = 0; d2[1] = 1;
  }
  else if (type == 'x') {
    d1[1] = 1;
    d2[2] = 1;
  }
  else if (type == 'y') {
    d1[2] = 1;
    d2[0] = 1;
  }
  else if (type == 'z') {
    d1[0] = 1;
    d2[1] = 1;
  }
  else {
    std::cout << "MxCartRect: invalid 'type' specification: " << type;
    throw 1;
  }

  this->edgeVertMatrix = MxDynMatrix<size_t>(this->numEdges, this->numVerts, 0);
  this->edgeVertMatrix(0, 0) = 1; this->edgeVertMatrix(0, 1) = 1;
  this->edgeVertMatrix(1, 0) = 1; this->edgeVertMatrix(1, 2) = 1;
  this->edgeVertMatrix(2, 3) = 1; this->edgeVertMatrix(2, 1) = 1;
  this->edgeVertMatrix(3, 3) = 1; this->edgeVertMatrix(3, 2) = 1;

  this->setLists();
}


template<size_t DIM>
std::vector<MxDimVector<double, DIM> > MxCartRect<DIM>::getVertices(const MxDimVector<double, DIM> & p) const {
  std::vector<MxDimVector<double, DIM> > res(4, p);

  res[0] += 0.5 * (-l1 * d1 - l2 * d2);
  res[1] += 0.5 * (l1 * d1 - l2 * d2);
  res[2] += 0.5 * (-l1 * d1 + l2 * d2);
  res[3] += 0.5 * (l1 * d1 + l2 * d2);

  return res;
}

#endif

