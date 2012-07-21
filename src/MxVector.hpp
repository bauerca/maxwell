#ifndef MX_VECTOR
#define MX_VECTOR

#include <vector>

#include "MxMultiVector.hpp"

#include "Epetra_Vector.h"
#include "Tpetra_Vector.hpp"

template<typename Scalar>
class MxVector : public MxMultiVector<Scalar> {
  public:
    MxVector(RCP<MxMap> map) : MxMultiVector<Scalar>(map, 1) {}

    MxVector(MxMultiVector<Scalar> const & mv, size_t vecIndex, bool deepcopy) :
        MxMultiVector<Scalar>(mv, std::vector<size_t>(1, vecIndex), deepcopy)
    {}

    /**
     *  Indices are always global here 
     */
    void replaceGlobalValue(MxIndex row, Scalar value) {
      MxMultiVector<Scalar>::replaceGlobalValue(row, 0, value);
    }

    void replaceLocalValue(MxIndex row, Scalar value) {
      MxMultiVector<Scalar>::replaceLocalValue(row, 0, value);
    }

    Scalar getValue(MxIndex localIndex) const {
      return MxMultiVector<Scalar>::operator()(0, localIndex);
    }

    RCP<Tpetra::Vector<Scalar, MxIndex> > getTpetraVector() const {
      return MxMultiVector<Scalar>::getTpetraMultiVector()->getVectorNonConst(0);
    }

    RCP<Epetra_Vector> getEpetraVector() const {
      return rcp(MxMultiVector<Scalar>::getEpetraMultiVector()->operator()(0));
    }

};

#endif
