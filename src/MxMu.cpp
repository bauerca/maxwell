#include "MxMu.hpp"


template<size_t DIM>
MxMu<DIM>::MxMu(std::string name) : mName(name),
mShape(0), mMu(MxDimMatrix<MxComplex, 3>::I()), mIsDiag(true) {}

template<size_t DIM>
MxMu<DIM>::MxMu(MxShape<DIM> const * shape, MxDimMatrix<MxComplex, 3> const & mu,
std::string name) : mName(name), mShape(shape), mMu(mu), mIsDiag(true) {
  mIsDiag = checkIsDiag(mMu);
} 

template<size_t DIM>
void MxMu<DIM>::setMu(MxShape<DIM> const * shape, MxDimMatrix<MxComplex, 3> const & mu) {
  mShape = shape;
  mMu = mu;
  mIsDiag = checkIsDiag(mMu);
}


template<size_t DIM>
bool MxMu<DIM>::checkIsDiag(MxDimMatrix<MxComplex, 3> const & mu) const {
  if (abs(mu(0, 1)) > MxUtil::dEps or
      abs(mu(0, 2)) > MxUtil::dEps or
      abs(mu(1, 2)) > MxUtil::dEps)
    return false;
  else
    return true;
}


template class MxMu<1>;
template class MxMu<2>;
template class MxMu<3>;
