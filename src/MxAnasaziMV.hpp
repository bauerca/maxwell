#ifndef MX_ANASAZI_MV
#define MX_ANASAZI_MV

#include "MxUtil.hpp"
#include "MxMap.hpp"
#include "MxMultiVector.hpp"

#include "AnasaziMultiVec.hpp"

template<typename Scalar>
class MxAnasaziMV : public Anasazi::MultiVec<Scalar>, public MxMultiVector<Scalar> {
  public:
    MxAnasaziMV(RCP<MxMap> map, size_t numVecs) : 
      MxMultiVector<Scalar>(map, numVecs) {}

    MxAnasaziMV(MxMultiVector<Scalar> const & mv) : 
      MxMultiVector<Scalar>(mv) {}

    MxAnasaziMV(MxMultiVector<Scalar> const & mv,
        std::vector<size_t> const & vecInds, bool deepcopy) : 
      MxMultiVector<Scalar>(mv, vecInds, deepcopy) {}

    virtual Anasazi::MultiVec<Scalar> * Clone(const int numVecs) const {
      //std::cout << "Clone called nvecs requested: " << numVecs << "\n";
      MxAnasaziMV<Scalar> * ptr = 
          new MxAnasaziMV<Scalar>(this->getMap(), numVecs);
      //std::cout << ptr->getRawMV()->NumVectors() << "\n";
      //std::cout << (*ptr->getRawMV())[0] << "\n";
      return ptr;
    }

    virtual Anasazi::MultiVec<Scalar> * CloneCopy() const {
      //std::cout << "CloneCopy called\n";
      MxAnasaziMV<Scalar> * ptr = 
          new MxAnasaziMV<Scalar>(*this);
      return ptr;
    }

    virtual Anasazi::MultiVec<Scalar> * CloneCopy(
        const std::vector<int> & index) const {
      std::vector<size_t> inds;
      MxUtil::convertVector(index, inds);
      //std::cout << "CloneCopy called with indices: ";
      //MxUtil::printStdVector(index);
      //MxUtil::printStdVector(inds);
      MxAnasaziMV<Scalar> * ptr = 
          new MxAnasaziMV<Scalar>(*this, inds, true); //deep copy
      return ptr;
    }

    virtual const Anasazi::MultiVec<Scalar> * CloneView(
        const std::vector<int> & index) const {
      std::vector<size_t> inds;
      MxUtil::convertVector(index, inds);
      //std::cout << "CloneView called with indices: ";
      //MxUtil::printStdVector(index);
      MxAnasaziMV<Scalar> * ptr = 
          new MxAnasaziMV<Scalar>(*this, inds, false); //shallow copy
      //std::cout << ptr->getRawMV()->NumVectors() << "\n";
      //std::cout << (*ptr->getRawMV())[0] << "\n";
      return ptr;
    }

    virtual Anasazi::MultiVec<Scalar> * CloneViewNonConst(
        const std::vector< int > &index) {
      std::vector<size_t> inds;
      MxUtil::convertVector(index, inds);
      //std::cout << "CloneViewNonConst called with indices: ";
      //MxUtil::printStdVector(index);
      MxAnasaziMV<Scalar> * ptr = 
          new MxAnasaziMV<Scalar>(*this, inds, false); //shallow copy
      //std::cout << ptr->getRawMV()->NumVectors() << "\n";
      //std::cout << (*ptr->getRawMV())[0] << "\n";
      return ptr;
    }

    virtual int GetVecLength() const {
      return int(this->getMap()->getGlobalNumIndices());
    }

    virtual int GetNumberVecs() const {
      return int(this->getNumVecs());
    }

    virtual void MvTimesMatAddMv(Scalar alpha,
        const Anasazi::MultiVec<Scalar> & A,
        const Teuchos::SerialDenseMatrix<int, Scalar> & B,
        Scalar beta);

    virtual void MvAddMv(Scalar alpha,
        const Anasazi::MultiVec<Scalar> & A, Scalar beta,
        const Anasazi::MultiVec<Scalar> & B);

    virtual void MvTransMv(Scalar alpha,
        const Anasazi::MultiVec<Scalar> & A,
        Teuchos::SerialDenseMatrix<int, Scalar> & B) const;

    virtual void MvDot(const Anasazi::MultiVec<Scalar> & A,
        std::vector<Scalar> & b) const {
      this->dot(dynamic_cast<const MxAnasaziMV<Scalar> &>(A), b);
    }

    virtual void MvNorm(std::vector<typename ScalarTraits<Scalar>::magnitudeType> & normvec) const {
      this->norm2(normvec);
    }

    virtual void SetBlock(const Anasazi::MultiVec<Scalar> & A,
        const std::vector<int> & index);

    virtual void MvScale(Scalar alpha) {
      //std::cout << "MvScale called\n";
      this->scale(alpha);
    }

    virtual void MvScale(const std::vector<Scalar> & alpha) {
      //std::cout << "MvScale called\n";
      this->scale(alpha);
    }

    virtual void MvRandom() {
      //std::cout << "MvRandom called\n";
      this->random();
    }

    virtual void MvInit(Scalar alpha) {
      //std::cout << "MvInit called\n";
      size_t n = this->getMap()->getNodeNumIndices();
      size_t nvecs = this->getNumVecs();
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < nvecs; ++j)
          this->replaceLocalValue(i, j, alpha);
    }

    virtual void MvPrint(std::ostream & os) const {
      os << "MxMultiVector\n";
    }

};

#endif
