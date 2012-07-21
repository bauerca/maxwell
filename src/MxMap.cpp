#include "MxMap.hpp"

#include "MxComm.hpp"
#include "Epetra_Map.h"

MxMap::MxMap(size_t numGlobalIndices, RCP<MxComm> comm) : 
mComm(comm) {
  setBaseMap(numGlobalIndices);
}

MxMap::MxMap(size_t numGlobalIndices, size_t numLocalIndices,
RCP<MxComm> comm) : 
mComm(comm) {
  setBaseMap(numGlobalIndices, numLocalIndices);
}

MxMap::MxMap(size_t numGlobalElements,
std::vector<MxIndex> const & globalElements,
RCP<MxComm> comm) : 
mComm(comm) {
  setBaseMap(numGlobalElements, globalElements);
}

MxMap::MxMap(std::vector<MxIndex> const & globalElements,
RCP<MxComm> comm) : 
mComm(comm) {
  setBaseMap(Teuchos::OrdinalTraits<size_t>::invalid(), globalElements);
}

#if 0
MxMap::MxMap(RCP<Tpetra::Map<MxIndex> const> map, RCP<MxComm> comm) : 
mComm(comm), mMap(map) {}
#endif

void MxMap::setBaseMap(size_t numGlobalIndices) {
#if LINALG_BASE == TPETRA
  mMap = rcp(new Tpetra::Map<MxIndex>(numGlobalIndices,
    0, mComm->getTeuchosComm()));
#elif LINALG_BASE == EPETRA
  mMap = rcp(new Epetra_Map(int(numGlobalIndices),
    0, *mComm->getEpetraComm()));
  setIntInds();
#endif
}

void MxMap::setBaseMap(size_t numGlobalIndices, size_t numLocalIndices) {
#if LINALG_BASE == TPETRA
  mMap = rcp(new Tpetra::Map<MxIndex>(numGlobalIndices, numLocalIndices,
    0, mComm->getTeuchosComm()));
#elif LINALG_BASE == EPETRA
  int numGlob = int(numGlobalIndices);
  if (numGlobalIndices == Teuchos::OrdinalTraits<size_t>::invalid())
    numGlob = -1;
  mMap = rcp(new Epetra_Map(numGlob, int(numLocalIndices), 0,
    *mComm->getEpetraComm()));
  setIntInds();
#endif
}

void MxMap::setBaseMap(size_t numGlobalIndices,
std::vector<MxIndex> const & globalIndices) {
#if LINALG_BASE == TPETRA
  mMap = rcp(new Tpetra::Map<MxIndex>(numGlobalIndices,
    Teuchos::ArrayView<MxIndex const>(globalIndices), 0,
    mComm->getTeuchosComm()));
#elif LINALG_BASE == EPETRA
  int numGlob = int(numGlobalIndices);
  if (numGlobalIndices == Teuchos::OrdinalTraits<size_t>::invalid())
    numGlob = -1;
  //int globalIndicesInt[globalIndices.size()];
  int * globalIndicesInt = new int[globalIndices.size()];
  for (size_t i = 0; i < globalIndices.size(); ++i)
    globalIndicesInt[i] = int(globalIndices[i]);
  mMap = rcp(new Epetra_Map(numGlob, int(globalIndices.size()),
    globalIndicesInt, 0, *mComm->getEpetraComm()));
  delete[] globalIndicesInt;
  setIntInds();
#endif
}


#if LINALG_BASE == EPETRA
void MxMap::setIntInds() {
  mGlobIndsInt = ArrayRCP<MxIndex>(mMap->NumMyElements());
  for (int i = 0; i < mMap->NumMyElements(); ++i)
    mGlobIndsInt[i] = MxIndex(mMap->MyGlobalElements()[i]);
}


RCP<Epetra_Map const> MxMap::getComplexEpetraMap() {
  if (mComplexEpetraMap == Teuchos::null) {
    int nGlob = 2*mMap->NumGlobalElements();
    int nLocal = 2*mMap->NumMyElements();

    int * inds = mMap->MyGlobalElements();
    int * indsC = new int[nLocal];
    for (int i = 0; i < mMap->NumMyElements(); ++i) {
      indsC[2*i] = 2*inds[i];
      indsC[2*i + 1] = 2*inds[i] + 1;
    }

    mComplexEpetraMap = rcp(new Epetra_Map(
      nGlob, nLocal, indsC, 0, *mComm->getEpetraComm()));
    delete[] indsC;
  }
  return mComplexEpetraMap;
}
#endif


int MxMap::getRemoteIndexList(std::vector<MxIndex> const & globInds,
std::vector<int> & nodeIDs, std::vector<MxIndex> & localInds) const {
  nodeIDs.resize(globInds.size());
  localInds.resize(globInds.size());

  return getRemoteIndexList(globInds.size(), &globInds[0],
    &nodeIDs[0], &localInds[0]);
}

int MxMap::getRemoteIndexList(size_t numInds, MxIndex const * globInds,
int * nodeIDs, MxIndex * localInds) const {

#if LINALG_BASE == TPETRA
  Tpetra::LookupStatus status;
  status = mMap->getRemoteIndexList(
    Teuchos::ArrayView<const MxIndex>(globInds, numInds),
    Teuchos::ArrayView<int>(nodeIDs, numInds),
    Teuchos::ArrayView<MxIndex>(localInds, numInds));

  if (status == Tpetra::AllIDsPresent)
    return 0;
  else
    return -1;
#elif LINALG_BASE == EPETRA
  int globIndsInt[numInds];
  int localIndsInt[numInds];
  for (size_t i = 0; i < numInds; ++i)
    globIndsInt[i] = int(globInds[i]);

  mMap->RemoteIDList(int(numInds), globIndsInt,
    nodeIDs, localIndsInt);

  for (size_t i = 0; i < numInds; ++i)
    localInds[i] = MxIndex(localIndsInt[i]);

  for (size_t i = 0; i < numInds; ++i)
    if (nodeIDs[i] == -1)
      return -1;

  return 0;

#endif
}


#if 0
template<typename Scalar>
RCP<Epetra_Map> MxMap::getEpetraMap() {
  if (mEpetraMap == Teuchos::null) {
    size_t numElems = mMap->getNodeNumElements();
    MxIndex const * elems = mMap->getNodeElementList().getRawPtr();

    if (Teuchos::ScalarTraits<Scalar>::isComplex) {
      std::vector<int> inds(2 * numElems);
      for (size_t j = 0; j < numElems; ++j) {
        inds[2*j] = int(2 * elems[j]);
        inds[2*j+1] = int(2 * elems[j] + 1);
      }

      mEpetraMap = rcp(new Epetra_Map(-1, int(inds.size()), &inds[0], 0,
          *mComm->getEpetraComm()));
    }
    else {
      std::vector<int> inds(numElems);
      for (size_t j = 0; j < numElems; ++j)
        inds[j] = int(elems[j]);

      mEpetraMap = rcp(new Epetra_Map(-1, int(numElems), &inds[0], 0,
          *mComm->getEpetraComm()));
    }
  }
  return mEpetraMap;
}

template RCP<Epetra_Map> MxMap::getEpetraMap<double>();
template RCP<Epetra_Map> MxMap::getEpetraMap<MxComplex>();
#endif
