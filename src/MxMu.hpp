#ifndef MX_MU
#define MX_MU

#include <vector>

#include "MxTypes.h"
#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxShape.hpp"

template<size_t DIM>
class MxMu {
  public:
  
    MxMu(std::string name);

    MxMu(MxShape<DIM> const * shape, MxDimMatrix<MxComplex, 3> const & mu,
      std::string name);

    std::string getName() const {return mName;}

    void setMu(MxShape<DIM> const * shape, MxDimMatrix<MxComplex, 3> const & mu);

    MxDimMatrix<MxComplex, 3> getMu() const {return mMu;}

    MxShape<DIM> const * getShape() const {return mShape;}

    bool isMuDiag() const {return mIsDiag;}


  private:
    std::string mName;

    bool mIsDiag;

    bool checkIsDiag(MxDimMatrix<MxComplex, 3> const & tensor) const;
    
    MxShape<DIM> const * mShape;

    MxDimMatrix<MxComplex, 3> mMu;

};


#endif
