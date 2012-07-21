#ifndef MX_MAP
#define MX_MAP

#include <vector>

#include "MxTypes.h"

#include "Tpetra_Map.hpp"
#include "Epetra_Map.h"

class MxComm;

template<typename Scalar> class MxCrsMatrix;
template<typename Scalar> class MxMultiVector;
template<typename Scalar> class MxVector;

class MxMap {
  public:
    /**
    Linear map creation
    */
    MxMap(size_t numGlobalIndices, RCP<MxComm> comm);

    /**
    Linear map creation with distribution control
    */
    MxMap(size_t numGlobalIndices, size_t numLocalIndices,
      RCP<MxComm> comm);

    MxMap(size_t numGlobalIndices,
      std::vector<MxIndex> const & myGlobalIndices, RCP<MxComm> comm);

    MxMap(std::vector<MxIndex> const & myGlobalIndices, RCP<MxComm> comm);

#if 0
    MxMap(RCP<Tpetra::Map<MxIndex> const> map, RCP<MxComm> comm);

    RCP<Tpetra::Map<MxIndex> const> getTpetraMap() const {
      return mMap;
    }

    /**
     *  FillComplete should have been called first
     */
    template<typename Scalar>
    RCP<Epetra_Map> getEpetraMap();
#endif

    size_t getNodeNumIndices() const {
#if LINALG_BASE == TPETRA
      return mMap->getNodeNumElements();
#elif LINALG_BASE == EPETRA
      return size_t(mMap->NumMyElements());
#endif
    }

    MxIndex const * getNodeIndexList() const {
#if LINALG_BASE == TPETRA
      return mMap->getNodeElementList().getRawPtr();
#elif LINALG_BASE == EPETRA
      return mGlobIndsInt.get();
#endif
    }

    MxIndex getGlobalIndex(MxIndex localIndex) const {
#if LINALG_BASE == TPETRA
      return mMap->getGlobalElement(localIndex);
#elif LINALG_BASE == EPETRA
      return mMap->GID(int(localIndex));
#endif
    }

    MxIndex getGlobalNumIndices() const {
#if LINALG_BASE == TPETRA
      return mMap->getGlobalNumElements();
#elif LINALG_BASE == EPETRA
      return mMap->NumGlobalElements();
#endif
    }

    MxIndex getLocalIndex(MxIndex globIndex) const {
#if LINALG_BASE == TPETRA
      return mMap->getLocalElement(globIndex);
#elif LINALG_BASE == EPETRA
      return mMap->LID(int(globIndex));
#endif
    }

    bool isNodeGlobalIndex(MxIndex globIndex) const {
#if LINALG_BASE == TPETRA
      return mMap->isNodeGlobalElement(globIndex);
#elif LINALG_BASE == EPETRA
      return mMap->MyGID(int(globIndex));
#endif
    }
    
    int getRemoteIndexList(std::vector<MxIndex> const & globInds,
        std::vector<int> & pids, std::vector<MxIndex> & localInds) const;

    int getRemoteIndexList(size_t numInds, MxIndex const * globInds,
        int * pids, MxIndex * localInds) const;

    RCP<MxComm> getComm() const {return mComm;}

    bool operator==(MxMap const & map) const {
      return map.mMap->SameAs(*mMap);
    }

    bool operator!=(MxMap const & map) const {
      return !map.mMap->SameAs(*mMap);
    }

    friend class MxCrsMatrix<double>;
    friend class MxMultiVector<double>;
    friend class MxVector<double>;
    friend class MxCrsMatrix<MxComplex>;
    friend class MxMultiVector<MxComplex>;
    friend class MxVector<MxComplex>;

  private:
    RCP<MxComm> mComm;

#if LINALG_BASE == TPETRA
    RCP<Tpetra::Map<MxIndex> const> mMap;
#elif LINALG_BASE == EPETRA
    RCP<Epetra_Map const> mMap, mComplexEpetraMap;

    ArrayRCP<MxIndex> mGlobIndsInt;

    void setIntInds();

    RCP<Epetra_Map const> getComplexEpetraMap();
#endif

    void setBaseMap(size_t numGlobalIndices);

    void setBaseMap(size_t numGlobalIndices,
      size_t numLocalIndices);

    void setBaseMap(size_t numGlobalIndices,
      std::vector<MxIndex> const & globalIndices);

};

#endif
