#ifndef MX_COMM
#define MX_COMM

#include "MxTypes.h"
#include "Tpetra_DefaultPlatform.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

class MxComm {
  public:
    MxComm() :
      mTeuchosComm(Tpetra::DefaultPlatform::getDefaultPlatform().getComm()),
#ifdef HAVE_MPI
      mEpetraComm(rcp(new Epetra_MpiComm(MPI_COMM_WORLD))) {}
#else
      mEpetraComm(rcp(new Epetra_SerialComm())) {}
#endif

    RCP<Epetra_Comm> getEpetraComm() const {return mEpetraComm;}

    RCP<Teuchos::Comm<int> const> getTeuchosComm() const {return mTeuchosComm;}

    int myPID() const {return mEpetraComm->MyPID();}

    size_t numProc() const {return mEpetraComm->NumProc();}

  private:
    RCP<Epetra_Comm> mEpetraComm;

    RCP<Teuchos::Comm<int> const> mTeuchosComm;

};

#endif
