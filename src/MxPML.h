#ifndef MX_PML
#define MX_PML

#include <complex>

#include "MxTypes.h"
#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxShape.hpp"
#include "MxGridField.hpp"

template<size_t DIM> class MxEMSim;

/**
 * Perfectly matched layer class. The PML implemented here is not equivalent
 * to the PML common to time-domain codes. To preserve the linearity of the
 * eigenvalue problem, this PML uses a conductivity that is linear in
 * frequency. In the PML region, the permittivity looks like
 *
 *   @f[ \varepsilon^{-1} =
 *      \frac{1}{\varepsilon_0 - i \sigma({\bf x},\omega)/\omega} @f]
 *
 * In time-domain codes, sigma is usually set to a constant value, giving
 * equivalent damping to all frequencies. The resulting frequency-dependence
 * of the permittivity complicates an eigensolve, so @f$ \sigma(\omega) @f$
 * is set to
 *
 *   @f[ \sigma(\omega) = \alpha({\bf x}) \omega @f]
 *
 */

template<size_t DIM>
class MxPML {
  public:
    MxPML(std::string name);

    void setShape(MxShape<DIM> const * shape) {mShape = shape;}

    /**
     * Specify the bounds of the conductivity profile. @f$ \alpha = 0 @f$
     * at start and @f$ \alpha_{\rm max} @f$ at end; varying according to
     * a power law (set with a call to setExponent(); default exponent is 2).
     *
     */
    void setProfileEndpoints(MxDimVector<double, DIM> start,
        MxDimVector<double, DIM> end);

    void setExponent(double exp) {mExp = exp;}

    void setAlphaMax(double alpha) {mAlphaMax = alpha;}

    double getAlpha(MxDimVector<double, DIM> point) const;

    MxDimMatrix<std::complex<double>, 3> getEps(MxDimVector<double, DIM> point) const;

    MxDimMatrix<std::complex<double>, 3> getInvEps(MxDimVector<double, DIM> point) const;

    MxDimMatrix<MxComplex, 3> getInvS(MxDimVector<double, DIM> point) const;

    MxDimMatrix<MxComplex, 3> getS(MxDimVector<double, DIM> point) const;

    const MxShape<DIM> * getShape() const {return mShape;}

    std::string getName() const {return mName;}
    
    bool isDiag() const {return mIsDiag;}

    static void alphaMap(MxEMSim<DIM> const * sim, MxGridField<DIM> const * field);

  private:
    std::string mName;

    double mAlphaMax;

    double mExp;

    double mWidth;

    bool mIsDiag;

    MxDimVector<double, DIM> mStart;
    MxDimVector<double, DIM> mEnd;
    MxDimVector<double, DIM> mDir;

    const MxShape<DIM> * mShape;

    // complex outer product of mDir
    MxDimMatrix<std::complex<double>, 3> mDirCOuter, mDirCOuterComplement;

    void setOuterProds();

    void setIsDiag();
};

template<size_t DIM>
MxPML<DIM>::MxPML(std::string name) : mName(name), mAlphaMax(0), mExp(2),
mWidth(0), mIsDiag(true), mStart(0), mEnd(0), mDir(0), mShape(0) {}


template<size_t DIM>
void MxPML<DIM>::setOuterProds() {
  MxDimVector<MxComplex, 3> mDirC(0);
  for (size_t i = 0; i < DIM; ++i)
    mDirC[i] = MxComplex(mDir[i], 0.0);

  mDirCOuter = MxDimMatrix<MxComplex, 3>(mDirC, mDirC);
  mDirCOuterComplement = MxDimMatrix<MxComplex, 3>::I() - mDirCOuter;
}

template<size_t DIM>
void MxPML<DIM>::setIsDiag() {
  int zeros = 0;
  for (size_t i = 0; i < DIM; ++i)
    if (abs(mDir[i]) < MxUtil::dEps)
      zeros++;

  if ((DIM == 3 and zeros == 2) or (DIM == 2 and zeros == 1))
    mIsDiag = true;
  else
    mIsDiag = false;
}

template<size_t DIM>
inline
void MxPML<DIM>::setProfileEndpoints(MxDimVector<double, DIM> start,
MxDimVector<double, DIM> end) {
  mStart = start;
  mEnd = end;
  mWidth = (end - start).norm();
  mDir = (end - start) / mWidth;
  setOuterProds();
  setIsDiag();
}


template<size_t DIM>
inline
MxDimMatrix<std::complex<double>, 3> MxPML<DIM>::getEps(MxDimVector<double, DIM> point) const {
  std::complex<double> s(1.0, mAlphaMax * pow((point - mStart).dot(mDir) / mWidth, mExp)); 
  return mDirCOuter / s + mDirCOuterComplement * s;
}

template<size_t DIM>
inline
MxDimMatrix<std::complex<double>, 3> MxPML<DIM>::getInvEps(MxDimVector<double, DIM> point) const {
  std::complex<double> s(1.0, mAlphaMax * pow((point - mStart).dot(mDir) / mWidth, mExp)); 
  return mDirCOuter * s + mDirCOuterComplement / s;
}

template<size_t DIM>
inline
MxDimMatrix<MxComplex, 3> MxPML<DIM>::getInvS(MxDimVector<double, DIM> point) const {
  MxComplex s(1.0, mAlphaMax * pow((point - mStart).dot(mDir) / mWidth, mExp)); 
  return mDirCOuter * s + mDirCOuterComplement;
}

template<size_t DIM>
inline
MxDimMatrix<MxComplex, 3> MxPML<DIM>::getS(MxDimVector<double, DIM> point) const {
  //std::cout << "pml width: " << mWidth << "\n";
  //std::cout << "pml exp: " << mExp << "\n";
  //std::cout << "point: " << point;
  //std::cout << "mstart: " << mStart;
  //std::cout << "alpha max: " << mAlphaMax << "\n";
  MxComplex s(1.0, mAlphaMax * pow((point - mStart).dot(mDir) / mWidth, mExp)); 
  return mDirCOuter / s + mDirCOuterComplement;
}

template<size_t DIM>
inline
double MxPML<DIM>::getAlpha(MxDimVector<double, DIM> point) const {
  return mAlphaMax * pow((point - mStart).dot(mDir) / mWidth, mExp);
}

#endif
