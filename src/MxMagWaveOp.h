#ifndef MX_MAG_WAVE_OP
#define MX_MAG_WAVE_OP

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <complex>

#include "MxEMSim.h"
#include "MxEMSimHierarchy.h"
#include "MxEMOps.h"

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Komplex_LinearProblem.h"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"
#include "AnasaziOperator.hpp"

template<size_t DIM, typename Scalar>
class MxMagWaveOp :
//public Epetra_Operator {
public Anasazi::Operator<Scalar> {

  public:
    MxMagWaveOp(RCP<MxEMSim<DIM> > sim);

    virtual ~MxMagWaveOp();

    //void setOperators(const MxEMSim<DIM> & aSim, MatMap & aMatMap);


    void eigValsToFreqs(std::vector<std::complex<double> > const & eigVals,
        std::vector<std::complex<double> > & freqs) const;

    MxMap const & getEMap() const {return *mEMap;}

    MxMap const & getBMap() const {return *mBMap;}

    MxMap const & getPsiMap() const {return *mPsiMap;}

    void magToElec(MxMultiVector<Scalar> const & mag,
        MxMultiVector<Scalar> & elec) const;

    void checkDivergences(MxMultiVector<Scalar> const & magVecs,
        std::vector<double> & residuals) const;

    void checkEigensolution(std::vector<std::complex<double> > const & eigVals,
        MxMultiVector<Scalar> const & eigVecs,
        std::vector<double> & residuals) const;

    void setShift(Scalar shift);

// Epetra_Operator interface
#if 0
    virtual int Apply(const Epetra_MultiVector & x, Epetra_MultiVector & y) const;

    virtual int ApplyInverse(const Epetra_MultiVector & x, Epetra_MultiVector & y) const;

    virtual double NormInf() const;

    virtual const Epetra_Comm & Comm() const;

    virtual const Epetra_Map & OperatorDomainMap() const;

    virtual const Epetra_Map & OperatorRangeMap() const;

    virtual int SetUseTranspose(bool use) {return 1;}

    virtual const char* Label() const {return "MxMagWaveOp";}

    virtual bool UseTranspose() const {return false;}

    virtual bool HasNormInf() const {return true;}
// Anasazi::Operator interface
#else
    virtual void Apply(Anasazi::MultiVec<Scalar> const & x,
        Anasazi::MultiVec<Scalar> & y) const;
#endif

  private:

    RCP<MxEMSim<DIM> > mSim;

    RCP<MxMap> mBMap, mEMap, mPsiMap;

    Teuchos::ParameterList mPList;

    int pid;

    void initMatrices();

    void initWorkVecs();

    void setLinearSolvers();

    bool mHasCurlNull;

    Epetra_CrsMatrix * getDiagBlocks(Epetra_CrsMatrix const & matrix) const;

    MxPolType pol;

    bool invert;

    int blockSize;

    double linTol;

    double shift;

    int linBasis;

    bool mIsComplex;

    bool mNeedEps;

    bool mNeedMu;

    mutable int mNumVecLinIters;
    mutable int mNumScaLinIters;
    mutable int mNumApplies;

    mutable double mVecCumSolveTime;
    mutable double mScaCumSolveTime;

    Teuchos::RCP<MxEMSimHierarchy<DIM> > simHier;

    RCP<AztecOO> mVecLaplSolver;

    RCP<AztecOO> mScaLaplSolver;

    RCP<Epetra_Operator> mVecLaplPrec, mScaLaplPrec;

    RCP<MxCrsMatrix<Scalar> > mVecLapl, mShiftedVecLapl, mPrecVecLapl;

    RCP<MxCrsMatrix<Scalar> > mScaLapl;

    RCP<MxCrsMatrix<Scalar> > mRhs;

    RCP<MxEMOps<DIM, Scalar> > mEMOps;

    RCP<Epetra_LinearProblem> mVecLaplProb, mScaLaplProb;

    mutable RCP<MxMultiVector<Scalar> > bWork1, bWork2;
    mutable RCP<MxMultiVector<Scalar> > eWork1;
    mutable RCP<MxMultiVector<Scalar> > psiWork1, psiWork2;

    RCP<MxMultiVector<double> > mVecCoords, mScaCoords;

    RCP<MxMultiVector<Scalar> > mVecNullspace, mScaNullspace;

    RCP<MxCrsMatrix<Scalar> > mNodeMatrix, mSMMatrix;

    void geoMGSetup();

    void amgSetup();

    void amgSetup(RCP<MxCrsMatrix<Scalar> > const & matrix,
        MxGridField<DIM> const & field,
        RCP<AztecOO> const & solver,
        RCP<Epetra_Operator> & prec,
        RCP<MxMultiVector<Scalar> > & nullspace,
        RCP<MxMultiVector<double> > & coords);

    void amgSetupRefMaxwell(RCP<MxCrsMatrix<Scalar> > const & matrix,
        MxGridField<DIM> const & field,
        RCP<AztecOO> const & solver,
        RCP<Epetra_Operator> & prec,
        RCP<MxMultiVector<Scalar> > & nullspace,
        RCP<MxMultiVector<double> > & coords);

    void iluSetup();

    void conditionDeyMittraVolumes(Epetra_CrsMatrix & dmInvVols, const Epetra_CrsMatrix & dmAreas, const MxEMSim<DIM> & theSim) const;

};


#endif
