#include "MxDielectric.hpp"


template<size_t DIM>
MxDielectric<DIM>::MxDielectric(std::string name) : mName(name),
mShape(0), mEps(MxDimMatrix<MxComplex, 3>::I()), mIsDiag(true) {}

template<size_t DIM>
MxDielectric<DIM>::MxDielectric(MxShape<DIM> const * shape, MxDimMatrix<MxComplex, 3> const & eps,
std::string name) : mName(name), mShape(shape), mEps(eps), mIsDiag(true) {
  if (!checkIsDiag(mEps))
    mIsDiag = false;
} 

template<size_t DIM>
void MxDielectric<DIM>::setEps(MxShape<DIM> const * shape, MxDimMatrix<MxComplex, 3> const & eps) {
  mShape = shape;
  mEps = eps;
  if (!checkIsDiag(mEps))
    mIsDiag = false;
}


template<size_t DIM>
bool MxDielectric<DIM>::checkIsDiag(MxDimMatrix<MxComplex, 3> const & eps) const {
  if (abs(eps(0, 1)) > MxUtil::dEps or
      abs(eps(0, 2)) > MxUtil::dEps or
      abs(eps(1, 2)) > MxUtil::dEps)
    return false;
  else
    return true;
}



template class MxDielectric<1>;
template class MxDielectric<2>;
template class MxDielectric<3>;
