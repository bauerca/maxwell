#ifndef MX_WAKE_SOLVER
#define MX_WAKE_SOLVER

#include "MxEMSim.h"
#include "MxPointCloud.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"

class Epetra_Operator;
class Epetra_MultiVector;

template<size_t DIM>
class MxWakeSolver {
  public:
    MxWakeSolver(const MxEMSim<DIM> * aSim, Teuchos::ParameterList pList);

    void solve();

    Teuchos::RCP<Epetra_MultiVector> getElecEigVecs(bool realPart = true) const;

    Teuchos::RCP<Epetra_MultiVector> getMagEigVecs(bool realPart = true) const;

    void saveFieldValsAtPoints(std::vector<MxDimVector<double, DIM> > const & points) const;

    void saveFieldValsAtPoints(MxPointCloud<DIM> const & pointCloud) const;

  private:
    typedef Epetra_Operator OP;
    typedef Epetra_MultiVector MV;

    int myPID;

    const MxEMSim<DIM> * sim;

    Teuchos::ParameterList pList, anasaziPList;

    Teuchos::RCP<OP> op;

    Teuchos::RCP<MV> initVec, uniFields;

    std::vector<std::complex<double> > eigVals, eigFreqs;

    std::vector<double> eigErrs;

    Teuchos::RCP<MV> eigVecs, elecEvecs;

    Teuchos::RCP<MV> magEvecsRe, magEvecsIm, elecEvecsRe, elecEvecsIm;

    Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > eigProb;

    Teuchos::RCP<Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> > eigMgr;

};

#endif
