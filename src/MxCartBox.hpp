#ifndef MX_CART_BOX
#define MX_CART_BOX

#include "MxDimVector.hpp"
#include "MxConvexPolyhedron.h"


class MxCartBox : public MxConvexPolyhedron {
  public:
    MxCartBox(double lx, double ly, double lz);
    MxCartBox(const MxDimVector<double, 3> & l);

    virtual ~MxCartBox() {}

    // don't include if you want to use the inherited method
    //virtual double volumeFraction(const MxShape<DIM> & aShape, const MxDimVector<double, DIM> & midPt) const;

    virtual double volume() const {return lx * ly * lz;}

    //virtual MxPolytope<DIM> intersection(const MxShape<DIM> & aShape) const;
    
    virtual std::vector<MxDimVector<double, 3> > getVertices(const MxDimVector<double, 3> & p) const;

  private:

    double lx, ly, lz;

    void setConnMatrices();

};


inline
std::vector<MxDimVector<double, 3> > MxCartBox::getVertices(const MxDimVector<double, 3> & p) const {
  std::vector<MxDimVector<double, 3> > res(8, p);
  
  res[0][0] -= 0.5 * lx; res[0][1] -= 0.5 * ly; res[0][2] -= 0.5 * lz;
  res[1][0] -= 0.5 * lx; res[1][1] -= 0.5 * ly; res[1][2] += 0.5 * lz;
  res[2][0] -= 0.5 * lx; res[2][1] += 0.5 * ly; res[2][2] -= 0.5 * lz;
  res[3][0] -= 0.5 * lx; res[3][1] += 0.5 * ly; res[3][2] += 0.5 * lz;
  res[4][0] += 0.5 * lx; res[4][1] -= 0.5 * ly; res[4][2] -= 0.5 * lz;
  res[5][0] += 0.5 * lx; res[5][1] -= 0.5 * ly; res[5][2] += 0.5 * lz;
  res[6][0] += 0.5 * lx; res[6][1] += 0.5 * ly; res[6][2] -= 0.5 * lz;
  res[7][0] += 0.5 * lx; res[7][1] += 0.5 * ly; res[7][2] += 0.5 * lz;

  return res;
}

#endif

