///////////////////////////
//
//  This class represents multiple dielectric objects. The order
//  in which dielectric objects are added determines which
//  object will 

#ifndef MX_DIELECTRIC
#define MX_DIELECTRIC

#include <vector>

#include "MxTypes.h"
#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxShape.hpp"

template<size_t DIM>
class MxDielectric {
  public:
  
    MxDielectric(std::string name);

    MxDielectric(MxShape<DIM> const * shape, MxDimMatrix<MxComplex, 3> const & eps,
      std::string name);

    std::string getName() const {return mName;}

    void setEps(MxShape<DIM> const * shape, MxDimMatrix<MxComplex, 3> const & eps);

    MxDimMatrix<MxComplex, 3> getEps() const {return mEps;}

    MxShape<DIM> const * getShape() const {return mShape;}

    bool isEpsDiag() const {return mIsDiag;}


  private:
    std::string mName;

    bool mIsDiag;

    bool checkIsDiag(MxDimMatrix<MxComplex, 3> const & tensor) const;
    
    MxShape<DIM> const * mShape;

    MxDimMatrix<MxComplex, 3> mEps;

};


#endif
