#ifndef MX_SOLVER
#define MX_SOLVER

#include "MxMultiVector.hpp"
#include "MxEMSim.h"
#include "MxPointCloud.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"

class Epetra_Operator;
class Epetra_MultiVector;

template<size_t DIM, typename Scalar>
class MxSolver {
  public:
    MxSolver(RCP<MxEMSim<DIM> > aSim, Teuchos::ParameterList pList);

    void solve();

#if 0
    Teuchos::RCP<Epetra_MultiVector> getElecEigVecs(bool realPart = true) const;

    Teuchos::RCP<Epetra_MultiVector> getMagEigVecs(bool realPart = true) const;

    void saveFieldValsAtPoints(std::vector<MxDimVector<double, DIM> > const & points) const;

    void saveFieldValsAtPoints(MxPointCloud<DIM> const & pointCloud) const;
#endif

  private:
    //typedef Epetra_Operator OP;
    //typedef Epetra_MultiVector MV;
    typedef Anasazi::Operator<Scalar> OP;
    typedef Anasazi::MultiVec<Scalar> MV;

    int myPID;

    RCP<MxEMSim<DIM> > sim;

    Teuchos::ParameterList pList, anasaziPList;

    Teuchos::RCP<OP> op;

    Teuchos::RCP<MV> mInitVec, uniFields;

    std::vector<MxComplex> mEigVals, mEigFreqs;

    std::vector<double> eigErrs;

    RCP<MxMultiVector<Scalar> > mBEigVecs, mEEigVecs;

    RCP<Anasazi::BasicEigenproblem<Scalar, MV, OP> > mEigProb;
    //RCP<Anasazi::BasicEigenproblem<double, MV, OP> > mEigProb; // for epetra interface

    RCP<Anasazi::SolverManager<Scalar, MV, OP> > mEigMgr;
    //RCP<Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> > mEigMgr; // for epetra

};

#endif
