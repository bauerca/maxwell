
#include "MxMagWaveOp.h"

#include <cmath>

#include "MxUtil.hpp"
#include "MxGridFieldIter.hpp"
#include "MxYeeDeyMittraCurlE.h"
#include "MxYeeDeyMittraCurlB.h"
#include "MxYeeDeyMittraDivB.h"
#include "MxYeeDeyMittraGradPsi.h"
#include "MxYeeDeyMittraFracs.h"
#include "MxYeeFitInvEps.h"
#include "MxYeeFitMu.h"
#include "MxYeeMitInvEps.h"
#include "MxEikxOp.h"
#include "MxAnasaziMV.hpp"
//#include "MxGeoMultigridPrec.h"

#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

#include "ml_MultiLevelPreconditioner.h"

template<size_t DIM, typename Scalar>
MxMagWaveOp<DIM, Scalar>::MxMagWaveOp(RCP<MxEMSim<DIM> > sim) : 
mSim(sim),
mBMap(sim->getField("bfield").getMap()),
mEMap(sim->getField("efield").getMap()),
pid(sim->getGrid().getComm()->myPID()),
mPList(sim->getParameters()),
mNeedEps(sim->hasDielectric() or sim->hasPML()),
mNeedMu(sim->hasMu() or sim->hasPML()),
mNumVecLinIters(0),
mNumScaLinIters(0),
mNumApplies(0) {
  Teuchos::ParameterList pList = mSim->getParameters();
  std::string metalAlg = pList.get("metal algorithm", "Dey-Mittra");
  pol = pList.get("polarization", TE);

  mIsComplex = pList.get("is complex", false);
  mHasCurlNull = (DIM == 3 or (DIM == 2 and pol == TM));
  if (mHasCurlNull)
    mPsiMap = sim->getField("psifield").getMap();

  invert = pList.get("eigensolver : invert", true);
  blockSize = pList.get("eigensolver : block size", 1);
  if (mIsComplex)
    blockSize *= 2;
  linTol = pList.get("eigensolver : tol", 1.0e-6);

  shift = pList.get("eigensolver : shift", 0.0);
  if (pid == 0) std::cout << "shift is: " << shift << "\n";


  initMatrices();

  // also sets linear solvers...
  setShift(shift);

  initWorkVecs();

  //saveOps(levelMats[0], "");

#if 0
  std::cout << "Operators in mag wave op\n";
  std::map<std::string, Teuchos::RCP<Epetra_CrsMatrix> >::const_iterator iter;
  for (iter = levelMats[0].begin(); iter != levelMats[0].end(); iter++) {
    std::cout << "  name: " <<
      iter->first << ", addr: " << iter->second.get() << "\n";
  }
#endif

  //print();
}


template<size_t DIM, typename Scalar>
MxMagWaveOp<DIM, Scalar>::~MxMagWaveOp() {
  if (pid == 0) {
    std::string invType = mPList.get("linear solver : type", "gmres");
    double vecIts = double(mNumVecLinIters)/double(mNumApplies);
    double scaIts = double(mNumScaLinIters)/double(mNumApplies);
    std::cout << "\n-----------MxMagWaveOp----------\n"
              << "  Total applications: " << mNumApplies << "\n"
              << "  VecLapl:\n"
              << "    - " << invType << " iterations: " << mNumVecLinIters << "\n"
              << "    - " << invType << " iters/app: " << vecIts << "\n"
              << "  ScaLapl:\n"
              << "    - " << invType << " iterations: " << mNumScaLinIters << "\n"
              << "    - " << invType << " iters/app: " << scaIts << "\n"
              << "--------------------------------\n\n";
  }

  // must destroy ml preconditioners first
  mScaLaplPrec = Teuchos::null;
  mVecLaplPrec = Teuchos::null;
}

#if 0
template<size_t DIM, typename Scalar>
void MxMagWaveOp<DIM, Scalar>::print() const {
  size_t numLevels = levelMats.size();
  MatMap::const_iterator matMapIter;
  for (size_t i = 0; i < numLevels; ++i) {
    if (pid == 0) std::cout << "Level " << i << " matrices: ";
    for (matMapIter = levelMats[i].begin(); matMapIter != levelMats[i].end(); ++matMapIter) {
      if (pid == 0) std::cout << matMapIter->first << ", ";
    }
    if (pid == 0) std::cout << "\n";
  }
}
#endif

template<size_t DIM, typename Scalar>
void MxMagWaveOp<DIM, Scalar>::initMatrices() {

  mEMOps = rcp(new MxEMOps<DIM, Scalar>(mSim));

  std::vector<std::string> names;

  // curl curl first
  names.push_back("curlB");
  if (mSim->hasDielectric())
    names.push_back("invEps");
  if (mSim->hasPEC())
    names.push_back("dmL");
  names.push_back("curlE");
  // fill complete it and store it
  mEMOps->multiplyOps(names, std::vector<bool>(names.size(), false),
    true, true, "curlCurl");

  // grad div if necessary
  if (mHasCurlNull) {
    names.clear();
    if (mSim->hasMu())
      names.push_back("mu");
    if (mSim->hasPEC())
      names.push_back("dmA");
    names.push_back("divB");
    if (mSim->hasPEC())
      names.push_back("dmVInv");
    if (mSim->hasMu())
      names.push_back("invMuVolAve");
    if (mSim->hasDielectric())
      names.push_back("invEpsVolAve");
    if (mSim->hasMu())
      names.push_back("invMuVolAve");
    names.push_back("gradPsi");
    if (mSim->hasMu())
      names.push_back("mu");
    if (mSim->hasPEC())
      names.push_back("dmA");
    // fill complete it store it
    mEMOps->multiplyOps(names, std::vector<bool>(names.size(), false),
      true, true, "gradDiv");
  }

  // if there is a curl nullspace, subtract gradDiv from curlCurl
  // to make vector laplacian
  if (mHasCurlNull) {
    names.clear();
    names.push_back("curlCurl");
    names.push_back("gradDiv");
    std::vector<Scalar> scalars;
    scalars.push_back(ScalarTraits<Scalar>::one());
    scalars.push_back(-ScalarTraits<Scalar>::one());

    mEMOps->addOps(names, scalars, true, "vecLapl");
    mVecLapl = mEMOps->getOp("vecLapl");

    // fill complete it store it
    //mVecLapl = mEMOps->addOps(names, scalars, true, "vecLapl");
    // fill complete it and don't store it
    //mVecLapl = mEMOps->addOps(names, scalars, false, "vecLapl");
    // Strip zeros from vector laplacian? there are many... Have to get
    // the function below working first
    //MxUtil::Epetra::stripZeros(*mVecLapl->getRawMatrix());
  } 
  else
    mVecLapl = mEMOps->getOp("curlCurl");

  // build scalar laplacian if necessary
  if (mHasCurlNull) {
    names.clear();

    names.push_back("gradPsi");
    if (mSim->hasMu())
      names.push_back("mu");
    if (mSim->hasPEC())
      names.push_back("dmA");
    names.push_back("divB");

    mScaLapl = mEMOps->multiplyOps(names, std::vector<bool>(names.size(), false),
      true, true, "scaLapl");
      //true, false, "scaLapl");
  }


  // The rhs of the generalized linear problem
  if (mSim->hasPEC() and mSim->hasMu()) {
    names.clear();
    names.push_back("mu");
    names.push_back("dmA");
    mRhs = mEMOps->multiplyOps(names,
      std::vector<bool>(names.size(), false), true);
  }
  else if (mSim->hasPEC()) {
    //std::cout << "pulling dmA as mRhs\n";
    mRhs = mEMOps->getOp("dmA");
  }
  else if (mSim->hasMu())
    mRhs = mEMOps->getOp("mu");
  else
    mRhs = rcp(new MxCrsMatrix<Scalar>(mBMap, ScalarTraits<Scalar>::one()));

  //RCP<MxCrsMatrix<Scalar> > mRhs2 = mRhs;
  //std::cout << mRhs2.has_ownership() << "\n\n\n\n";
  mEMOps->saveOps();
}

template<size_t DIM, typename Scalar>
void MxMagWaveOp<DIM, Scalar>::setShift(Scalar shift) {
  
  //mShiftedVecLapl = rcp(new MxCrsMatrix<Scalar>(mBMap, mBMap));
  mShiftedVecLapl = rcp(new MxCrsMatrix<Scalar>(mBMap));

  MxCrsMatrix<Scalar>::matrixMatrixAdd(*mVecLapl, false,
    ScalarTraits<Scalar>::one(), *mRhs, false, -shift,
    *mShiftedVecLapl, true);
  mShiftedVecLapl->save("shiftedVecLapl");

  // need to also update the matrix used to precondition
  Scalar precShift;
  //precShift = (-ScalarTraits<Scalar>::one() + MxUtil::i<Scalar>()) * shift;
  if (ScalarTraits<Scalar>::isComplex)
    precShift = MxUtil::i<Scalar>() * shift;
  else
    precShift = shift;

  mPrecVecLapl = rcp(new MxCrsMatrix<Scalar>(mBMap));
  MxCrsMatrix<Scalar>::matrixMatrixAdd(*mVecLapl, false,
    ScalarTraits<Scalar>::one(), *mRhs, false, precShift,
    *mPrecVecLapl, true);
  mPrecVecLapl->save("precVecLapl");
  
  setLinearSolvers();
}

template<size_t DIM, typename Scalar>
void MxMagWaveOp<DIM, Scalar>::initWorkVecs() {
  bWork1 = rcp(new MxMultiVector<Scalar>(mBMap, blockSize));

  if (mHasCurlNull) {
    psiWork1 = rcp(new MxMultiVector<Scalar>(mPsiMap, blockSize));
    psiWork2 = rcp(new MxMultiVector<Scalar>(mPsiMap, blockSize));
  }
}

template<size_t DIM, typename Scalar>
void MxMagWaveOp<DIM, Scalar>::setLinearSolvers() {

  mVecLaplProb = rcp(new Epetra_LinearProblem(
      mShiftedVecLapl->getRawMatrix().get(), NULL, NULL));
  mVecLaplSolver = rcp(new AztecOO(*mVecLaplProb));

  if (mHasCurlNull) {
    mScaLaplProb = rcp(new Epetra_LinearProblem(
        mScaLapl->getRawMatrix().get(), NULL, NULL));
    mScaLaplSolver = rcp(new AztecOO(*mScaLaplProb));
  }


  // get the type of the preconditioner
  std::string precType = mPList.get("linear solver : prec type", "none");
  if (precType == "gmg") {
    if (pid == 0) std::cout << "Using GMRES with geometric multigrid preconditioner.\n";
    if (pid == 0) std::cout << "Not implemented.\n";
    exit(EXIT_FAILURE);
    //geoMGSetup();
  }
  else if (precType == "amg") {
    if (pid == 0) std::cout << "Using GMRES with algebraic multigrid preconditioner.\n";
    amgSetup();
  }
  else if (precType == "ilu") {
    if (pid == 0) std::cout << "Using GMRES with ILU preconditioner.\n";
    iluSetup();
  }
  else if (precType == "none") {
    if (pid == 0) std::cout << "MxMagWaveOp: WARNING: no preconditioner!\n";
  }
  else {
    std::cout << "MxMagWaveOp: Unrecognized preconditioner type '" << precType
      << "'.\n";
    exit(EXIT_FAILURE);
  }

  // get the type of the linear solver
  std::string invType = mPList.get("linear solver : type", "gmres");
  int linSolver;
  if (invType == "gmres")
    linSolver = AZ_gmres;
  else if (invType == "cg")
    linSolver = AZ_cg;
  else if (invType == "bicgstab")
    linSolver = AZ_bicgstab;
  else {
    std::cout << "MxMagWaveOp: solver type '" << invType << "' not implemented.\n";
    throw 1;
  }

  linBasis = mPList.get("linear solver : basis", 20);
  int linOutput = mPList.get("linear solver : verbosity", 1);

  mVecLaplSolver->SetAztecOption(AZ_output, linOutput);
  mVecLaplSolver->SetAztecOption(AZ_kspace, linBasis);
  mVecLaplSolver->SetAztecOption(AZ_solver, linSolver);
  mVecLaplSolver->SetAztecOption(AZ_conv, AZ_noscaled);

  if (mHasCurlNull) {
    mScaLaplSolver->SetAztecOption(AZ_output, linOutput);
    mScaLaplSolver->SetAztecOption(AZ_kspace, linBasis);
    mScaLaplSolver->SetAztecOption(AZ_solver, linSolver);
    mScaLaplSolver->SetAztecOption(AZ_conv, AZ_noscaled);
  }

}

template<size_t DIM, typename Scalar>
void MxMagWaveOp<DIM, Scalar>::iluSetup() {

  double dropTol = mPList.get("linear solver : ilu : drop tol", 1.e-3);
  double fill = mPList.get("linear solver : ilu : fill", 100.0);

  double condest;

  mVecLaplSolver->SetPrecMatrix(mPrecVecLapl->getRawMatrix().get());

  mVecLaplSolver->SetAztecOption(AZ_precond, AZ_dom_decomp);
  mVecLaplSolver->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
  mVecLaplSolver->SetAztecOption(AZ_keep_info, 1);
  mVecLaplSolver->SetAztecParam(AZ_drop, dropTol);
  mVecLaplSolver->SetAztecParam(AZ_ilut_fill, fill);
  mVecLaplSolver->ConstructPreconditioner(condest);
  mVecLaplSolver->SetAztecOption(AZ_pre_calc, AZ_reuse);

  if (mHasCurlNull) {
    mScaLaplSolver->SetAztecOption(AZ_precond, AZ_dom_decomp);
    mScaLaplSolver->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    mScaLaplSolver->SetAztecOption(AZ_keep_info, 1);
    mScaLaplSolver->SetAztecParam(AZ_drop, dropTol);
    mScaLaplSolver->SetAztecParam(AZ_ilut_fill, fill);
    mScaLaplSolver->ConstructPreconditioner(condest);
    mScaLaplSolver->SetAztecOption(AZ_pre_calc, AZ_reuse);
  }

}




template<size_t DIM, typename Scalar>
void MxMagWaveOp<DIM, Scalar>::amgSetup() {

  // vector laplacian
  amgSetup(mPrecVecLapl, mSim->getField("bfield"), mVecLaplSolver,
    mVecLaplPrec, mVecNullspace, mVecCoords);

  if (mHasCurlNull) {
    // scalar laplacian
    amgSetup(mScaLapl, mSim->getField("psifield"), mScaLaplSolver,
      mScaLaplPrec, mScaNullspace, mScaCoords);
  }
}

template<size_t DIM, typename Scalar>
void MxMagWaveOp<DIM, Scalar>::amgSetup(
RCP<MxCrsMatrix<Scalar> > const & matrix,
MxGridField<DIM> const & field,
RCP<AztecOO> const & solver,
RCP<Epetra_Operator> & prec,
RCP<MxMultiVector<Scalar> > & nullspace,
RCP<MxMultiVector<double> > & coords) {
  int amgSweeps = mPList.get("linear solver : amg : sweeps", 1);
  std::string amgSmootherType = mPList.get("linear solver : amg : smoother", "Chebyshev");
  std::string coarseSmootherType = mPList.get("linear solver : amg : coarse smoother", "Amesos-KLU");
  int maxLevels = mPList.get("linear solver : amg : max levels", 10);

  if (coarseSmootherType == "Amesos")
    coarseSmootherType = "Amesos-KLU";

  Teuchos::ParameterList mlList;
  ML_Epetra::SetDefaults("SA", mlList);

  mlList.set("ML output", 10);
  mlList.set("max levels", maxLevels);

  //mlList.set("coarse: max size", 30);
  //mlList.set("eigen-analysis: type", "power-method");
  //mlList.set("eigen-analysis: iterations", 100);
  mlList.set("cycle applications", 1);
  mlList.set("prec type", "full-MGV");
  
  mlList.set("aggregation: type", "Uncoupled");
  //mlList.set("aggregation: type", "MIS");
  
  //mlList.set("smoother: type", "Gauss-Seidel");
  mlList.set("smoother: type", amgSmootherType);
  mlList.set("smoother: sweeps", amgSweeps);
  //mlList.set("smoother: damping factor", 2./3.);
  
  mlList.set("coarse: type", coarseSmootherType);
  mlList.set("coarse: sweeps", amgSweeps);

  //mlList.set("energy minimization: enable", true);

  nullspace = MxGridField<DIM>::template uniformFields<Scalar>(field);

  // should test whether applying eikx operator here helps.
  // for now leave the nullspace as the constant vectors

  // another possible improvement here... untested
#if 0
  // zero part of uniform vector in pml
  MxGridFieldIter<DIM> iter(&mSim->getField("bfield"));
  int row;
  bool skip;
  for (iter.begin(); !iter.atEnd(); iter.bump()) {
    row = iter.getGlobCompIndx();
    skip = false;
    for (size_t i = 0; i < mSim->numPMLs(); ++i) {
      if (mSim->getPML(i)->getShape()->func(iter.getCoord()) > 0)
        skip = true;
    }

    for (int i = 0; i < nVec; ++i)
      uni.ReplaceGlobalValue(row, i, 0.0);
  }
#endif


  //std::cout << *vecNullspace;
  double * nullspaceView;
  int jump;
  nullspace->getRawMV()->ExtractView(&nullspaceView, &jump);
  //std::cout << "nullspace vector pointer: " << vecNullspaceView << "\n";

  if (ScalarTraits<Scalar>::isComplex and nullspace->getNumVecs() > 3)
    mlList.set("PDE equations", 2);
  mlList.set("null space: type", "pre-computed");
  mlList.set("null space: vectors", nullspaceView);
  mlList.set("null space: dimension", int(nullspace->getNumVecs()));
  mlList.set("null space: add default vectors", false);

  // do matrix repartitioning with Zoltan
  mlList.set("repartition: enable", 1);
  mlList.set("repartition: Zoltan dimensions", int(DIM));

  // if Scalar is complex, this multivector will have half as many
  // entries as a field vector
  coords = field.coords();
  mlList.set("x-coordinates", (*coords->getRawMV())[0]);
  if (DIM > 1) mlList.set("y-coordinates", (*coords->getRawMV())[1]);
  if (DIM > 2) mlList.set("z-coordinates", (*coords->getRawMV())[2]);


  //mlList.set("max levels", 2);
  mlList.set("coarse: type", "Chebyshev");
  mlList.set("coarse: sweeps", 4);

  // forming multigrid preconditioner based on whole input matrix
  // There is some evidence that preconditioning on only the diagonal
  // blocks for complex-equivalent real matrices is better.
  prec = rcp(new ML_Epetra::MultiLevelPreconditioner(
    *matrix->getRawMatrix(), mlList));
  solver->SetPrecOperator(prec.get());

}

// input matrix is one created with Komplex. Returns the real part, or the diagonal
// blocks.
template<size_t DIM, typename Scalar>
Epetra_CrsMatrix * MxMagWaveOp<DIM, Scalar>::getDiagBlocks(Epetra_CrsMatrix const & matrix) const {
  Epetra_CrsMatrix * res = new Epetra_CrsMatrix(Copy, matrix.RowMap(), 0);

  int max = matrix.MaxNumEntries();
  double * vals;
  int * cols;
  int col;

  int numEntries;
  int mod2;

  int * globElems = matrix.RowMap().MyGlobalElements();
  int * globColElems = matrix.ColMap().MyGlobalElements();
  int numMyElems = matrix.RowMap().NumMyElements();
  for (int i = 0; i < numMyElems; ++i) {
    mod2 = globElems[i] % 2;
    matrix.ExtractMyRowView(i, numEntries, vals, cols);
    //std::cout << "row " << globElems[i] << "; ";

    for (int j = 0; j < numEntries; ++j) {
      col = globColElems[cols[j]];
      //std::cout << col << ", ";
      if ((col % 2) == mod2)
        res->InsertGlobalValues(globElems[i], 1, &vals[j], &col);
    }
    //std::cout << "\n";
  }

  res->FillComplete(matrix.DomainMap(), matrix.RangeMap());
  return res;
}

#if 0
template<size_t DIM>
void MxMagWaveOp<DIM>::geoMGSetup() {
  Teuchos::ParameterList pList = mSim->getParameters();

  int targetNumLevels = pList.get("linear solver : levels", 1);

  mSimHier = Teuchos::rcp(new MxEMSimHierarchy<DIM>(mSim, targetNumLevels));

  size_t numLevels = mSimHier->getNumLevels();

  levelMats.resize(numLevels);

  // inputs to MxGeoMultigridPrec
  std::vector<Teuchos::RCP<Epetra_CrsMatrix> > vecLapls, scaLapls;
  std::vector<Teuchos::RCP<Epetra_CrsMatrix> > vecRHSCoarseners, scaRHSCoarseners;
  std::vector<Teuchos::RCP<Epetra_CrsMatrix> > vecLHSCoarseners, scaLHSCoarseners;
  std::vector<Teuchos::RCP<Epetra_CrsMatrix> > vecRHSRefiners, scaRHSRefiners;
  std::vector<Teuchos::RCP<Epetra_CrsMatrix> > vecLHSRefiners, scaLHSRefiners;

  // get pointers to finest level laplacian matrices
  vecLapls.push_back(levelMats[0].find("vecLapl")->second);
  if (DIM == 3 or (DIM == 2 and pol == "TM"))
    scaLapls.push_back(levelMats[0].find("scaLapl")->second);

  // get pointers to coarse level laplacian matrices
  const MxEMSim<DIM> * coarseLevelSim, * fineLevelSim = mSim;
  for (size_t level = 1; level < numLevels; ++level) {
    if (pid == 0) std::cout << "Initializing coarse level " << level << "\n";
    coarseLevelSim = mSimHier->getSim(level);

    if (pid == 0) std::cout << "Setting operators\n";
    if (mIsComplex)
      setSingleOpsComplex(*coarseLevelSim, levelMats[level]);
    else
      setSingleOps(*coarseLevelSim, levelMats[level], false);
    combineOperators(*coarseLevelSim, levelMats[level]);

    const MxGridField<DIM> & coarseB = coarseLevelSim->getField("bfield");
    const MxGridField<DIM> & fineB = fineLevelSim->getField("bfield");

    if (pid == 0) std::cout << "Setting B-field coarsener...";
    vecLHSCoarseners.push_back(Teuchos::rcp(new Epetra_CrsMatrix(fineB.getInterpolator(coarseB))));
    vecRHSCoarseners.push_back(Teuchos::rcp(new Epetra_CrsMatrix(fineB.getInterpolator(coarseB))));
    if (pid == 0) std::cout << "Setting B-field refiner...";
    vecLHSRefiners.push_back(Teuchos::rcp(new Epetra_CrsMatrix(coarseB.getInterpolator(fineB))));
    vecRHSRefiners.push_back(Teuchos::rcp(new Epetra_CrsMatrix(coarseB.getInterpolator(fineB))));

    vecLapls.push_back(levelMats[level].find("vecLapl")->second);

    if (DIM == 3 or (DIM == 2 and pol == "TM")) {
   
      const MxGridField<DIM> & coarsePsi = coarseLevelSim->getField("psifield");
      const MxGridField<DIM> & finePsi = fineLevelSim->getField("psifield");

      if (pid == 0) std::cout << "Setting Psi-field coarsener...";
      scaRHSCoarseners.push_back(Teuchos::rcp(new Epetra_CrsMatrix(finePsi.getInterpolator(coarsePsi))));
      if (pid == 0) std::cout << "done.\n";

      if (pid == 0) std::cout << "Setting Psi-field refiner...";
      scaRHSRefiners.push_back(Teuchos::rcp(new Epetra_CrsMatrix(coarsePsi.getInterpolator(finePsi))));
      if (pid == 0) std::cout << "done.\n";

      scaLapls.push_back(levelMats[level].find("scaLapl")->second);
    }

    fineLevelSim = coarseLevelSim;
  }

  MxGeoMultigridPrec * vecMGPtr = new MxGeoMultigridPrec(pList);

  vecMGPtr->setOperators(vecLapls);
  vecMGPtr->setRHSCoarseners(vecRHSCoarseners);
  vecMGPtr->setLHSCoarseners(vecLHSCoarseners);
  vecMGPtr->setRHSRefiners(vecRHSRefiners);
  vecMGPtr->setLHSRefiners(vecLHSRefiners);
  vecMGPtr->setup();
  vecLaplPrec = Teuchos::rcp(vecMGPtr);
  mVecLaplSolver->SetPrecOperator(vecLaplPrec.get());

  if (DIM == 3 or (DIM == 2 and pol == "TM")) {
    Teuchos::ParameterList scaLaplList(pList);
    scaLaplList.set("geo-mg : remove const field", true);
    MxGeoMultigridPrec * scaMGPtr = new MxGeoMultigridPrec(scaLaplList);

    scaMGPtr->setOperators(scaLapls);
    scaMGPtr->setCoarseners(scaRHSCoarseners);
    scaMGPtr->setRefiners(scaRHSRefiners);
    scaMGPtr->setup();
    scaLaplPrec = Teuchos::rcp(scaMGPtr);
    mScaLaplSolver->SetPrecOperator(scaLaplPrec.get());
  }

}


template<size_t DIM>
void MxMagWaveOp<DIM>::saveOps(const MatMap & ops, std::string suffix) const {
  if (pid == 0) std::cout << "Saving operators in MxMagWaveOp...\n";
  std::map<std::string, Teuchos::RCP<Epetra_CrsMatrix> >::const_iterator iter;
  std::string name;
  for (iter = ops.begin(); iter != ops.end(); iter++) {
    name = iter->first + suffix + std::string(".mm");
    EpetraExt::RowMatrixToMatrixMarketFile(name.c_str(), *iter->second);
  }
}
#endif


// Anasazi::Operator interface
#if 1
template<size_t DIM, typename Scalar>
void MxMagWaveOp<DIM, Scalar>::Apply(const Anasazi::MultiVec<Scalar> & x,
Anasazi::MultiVec<Scalar> & y) const {
  MxAnasaziMV<Scalar> const & x2 = dynamic_cast<MxAnasaziMV<Scalar> const &>(x);
  MxAnasaziMV<Scalar> & y2 = dynamic_cast<MxAnasaziMV<Scalar> &>(y);

  if (pid == 0) std::cout << "Applying MagWaveOp!\n";
  mNumApplies++;

  if (invert) {
    mRhs->apply(x2, y2);

    // vec lapl solver is epetra-based
    mVecLaplSolver->SetRHS(y2.getRawMV().get());
    bWork1->set(ScalarTraits<Scalar>::zero());
    mVecLaplSolver->SetLHS(bWork1->getRawMV().get());

    // do first inversion: Y = (M - shift*A)^{-1}*X
    double trueRes = 1;
    while (trueRes > linTol) {
      mVecLaplSolver->recursiveIterate(linBasis, linTol);
      if (pid == 0) std::cout << "Getting true residual...\n";
      trueRes = mVecLaplSolver->TrueResidual();
      mNumVecLinIters += mVecLaplSolver->NumIters();
      if (pid == 0) std::cout << "True residual of vector Laplacian solve: " << trueRes << "\n";
    }
    if (pid == 0) std::cout << "True residual of vector Laplacian solve: " << trueRes << "\n";

    if (mHasCurlNull) {
      // now do projection
      mRhs->apply(*bWork1, y2);
      mEMOps->getOp("divB")->apply(y2, *psiWork1);

      // remove constant psi field from psiWork1
      //scaLaplacianPrec->removeConstField(*psiWork1);
      //MxUtil::Epetra::removeConstField(*psiWork1);

      // sca lapl solver is epetra-based
      mScaLaplSolver->SetRHS(psiWork1->getRawMV().get());
      psiWork2->set(ScalarTraits<Scalar>::zero());
      mScaLaplSolver->SetLHS(psiWork2->getRawMV().get());

      trueRes = 1;
      while (trueRes > linTol) {
        mScaLaplSolver->recursiveIterate(linBasis, linTol);
        trueRes = mScaLaplSolver->TrueResidual();
        mNumScaLinIters += mScaLaplSolver->NumIters();
      }
      if (pid == 0) std::cout << "True residual of scalar Laplacian solve: " << trueRes << "\n";

      //levelMats[0].find("divB")->second->Multiply(true, *psiWork2, y);
      //levelMats[0].find("gradPsi")->second->Multiply(false, *psiWork2, y);
      mEMOps->getOp("gradPsi")->apply(*psiWork2, y2);

      y2.update(ScalarTraits<Scalar>::one(), *bWork1,
          -ScalarTraits<Scalar>::one());
    }
    else {
      y2 = *bWork1;
      if (pid == 0) std::cout << "Get here?\n";
    }
  }
  // else if no inversion
  else {
#if 0
    levelMats[0].find("curlB")->second->Multiply(false, x, *eWork1);
    if (mSim->hasPEC())
      levelMats[0].find("dmL")->second->Multiply(false, *eWork1, *eWork1);
    levelMats[0].find("curlE")->second->Multiply(false, *eWork1, y);
    if (mSim->hasPEC())
      levelMats[0].find("dmA")->second->ApplyInverse(y, y);
    y.Update(-shift, x, 1.0);
#endif
    if (pid == 0) std::cout << "Non inverse Mag wave op not implemented\n";
  }
}
#endif // Anasazi interface


// Epetra_Operator interface
#if 0
template<size_t DIM, typename Scalar>
int MxMagWaveOp<DIM, Scalar>::Apply(const Epetra_MultiVector & x, Epetra_MultiVector & y) const {
  if (invert) {

    mRhs->getRawMatrix()->Apply(x, y);

    mVecLaplSolver->SetRHS(&y);
    //mVecLaplProb->SetRHS(&y);
    bWork1->getRawMV()->PutScalar(0.0);
    mVecLaplSolver->SetLHS(bWork1->getRawMV().get());
    //mVecLaplProb->SetLHS(bWork1.get());

    // do first inversion: Y = (M - shift*A)^{-1}*X
    double trueRes = 1;
    while (trueRes > linTol) {
      mVecLaplSolver->recursiveIterate(linBasis, linTol);
      trueRes = mVecLaplSolver->TrueResidual();
      //num_lin_iters1_ += solver_->NumIters();
    }
    if (pid == 0) std::cout << "True residual of vector Laplacian solve: " << trueRes << "\n";

    if (mHasCurlNull) {
      // now do projection
      mRhs->getRawMatrix()->Apply(*bWork1->getRawMV(), y);
      mEMOps->getOp("divB")->getRawMatrix()->Apply(y, *psiWork1->getRawMV());

      // remove constant psi field from psiWork1
      //scaLaplacianPrec->removeConstField(*psiWork1);
      //MxUtil::Epetra::removeConstField(*psiWork1);

      mScaLaplSolver->SetRHS(psiWork1->getRawMV().get());
      psiWork2->getRawMV()->PutScalar(0.0);
      mScaLaplSolver->SetLHS(psiWork2->getRawMV().get());

      trueRes = 1;
      while (trueRes > linTol) {
        mScaLaplSolver->recursiveIterate(linBasis, linTol);
        trueRes = mScaLaplSolver->TrueResidual();
      }
      if (pid == 0) std::cout << "True residual of scalar Laplacian solve: " << trueRes << "\n";

      //levelMats[0].find("divB")->second->Multiply(true, *psiWork2, y);
      //levelMats[0].find("gradPsi")->second->Multiply(false, *psiWork2, y);
      mEMOps->getOp("gradPsi")->getRawMatrix()->Apply(*psiWork2->getRawMV(), y);

      y.Update(1.0, *bWork1->getRawMV(), -1.0);
    }
    else
      y = *bWork1->getRawMV();
  }
  // else if no inversion
  else {
#if 0
    levelMats[0].find("curlB")->second->Multiply(false, x, *eWork1);
    if (mSim->hasPEC())
      levelMats[0].find("dmL")->second->Multiply(false, *eWork1, *eWork1);
    levelMats[0].find("curlE")->second->Multiply(false, *eWork1, y);
    if (mSim->hasPEC())
      levelMats[0].find("dmA")->second->ApplyInverse(y, y);
    y.Update(-shift, x, 1.0);
#endif
    if (pid == 0) std::cout << "Non inverse Mag wave op not implemented\n";
  }

  return 0;
}

template<size_t DIM, typename Scalar>
int MxMagWaveOp<DIM, Scalar>::ApplyInverse(const Epetra_MultiVector & x, Epetra_MultiVector & y) const {
  if (pid == 0) std::cout << "MxMagWaveOp::ApplyInverse: does nothing.\n";
  return -1;
}

template<size_t DIM, typename Scalar>
double MxMagWaveOp<DIM, Scalar>::NormInf() const {
  return mEMOps->getOp("vecLapl")->getRawMatrix()->NormInf();
}

template<size_t DIM, typename Scalar>
const Epetra_Comm & MxMagWaveOp<DIM, Scalar>::Comm() const {
  return *mSim->getGrid().getComm()->getEpetraComm();
}

template<size_t DIM, typename Scalar>
const Epetra_Map & MxMagWaveOp<DIM, Scalar>::OperatorDomainMap() const {
  return mEMOps->getOp("vecLapl")->getRawMatrix()->OperatorDomainMap();
  //if (mIsComplex)
  //  return mVecLaplProbK->KomplexProblem()->GetOperator()->OperatorDomainMap();
  //else
  //  return mSim->getField("bfield").getMap();
}

template<size_t DIM, typename Scalar>
const Epetra_Map & MxMagWaveOp<DIM, Scalar>::OperatorRangeMap() const {
  return mEMOps->getOp("vecLapl")->getRawMatrix()->OperatorRangeMap();
  //if (mIsComplex)
  //  return mVecLaplProbK->KomplexProblem()->GetOperator()->OperatorRangeMap();
  //else
  //  return mSim->getField("bfield").getMap();
}

#endif // Epetra_Operator interface


#if 0
template<size_t DIM>
void MxMagWaveOp<DIM>::conditionDeyMittraVolumes(Epetra_CrsMatrix & dmInvVols, const Epetra_CrsMatrix & dmAreas, const MxEMSim<DIM> & theSim) const {
  const MxGridField<DIM> & psiField = theSim.getField("psifield");
  const MxGridField<DIM> & bField = theSim.getField("bfield");

  MxGridFieldIter<DIM> iter(&psiField);

  size_t numBComps = bField.getNumComps();
  int N = 2 * numBComps;

  MxDimVector<int, DIM> cell;
  double aprods[2 * numBComps];
  double a1, a2;
  double rowsum;
  double vfrac;
  int comp, bump;
  char comps[3] = {'x', 'y', 'z'};

  double large = 1e2;
  int largeRatios = 0;

  for (iter.begin(); !iter.atEnd(); iter.bump()) {
    cell = iter.getCell();

    vfrac = psiField.getCompFrac(0, cell, "pec");

    if (vfrac != 1.0) {

      std::cout << "Cell "; cell.print();
      std::cout << "  vfrac = " << vfrac << "\n";
      for (int i = 0; i < N; ++i) {
        comp = i / 2;
        bump = i % 2;
        cell[comp] += bump;
        a1 = bField.getCompFrac(comp, cell, "pec");
        cell[comp] -= bump;

        rowsum = 0;

        for (int j = 0; j < N; ++j) {
          comp = j / 2;
          bump = j % 2;
          cell[comp] += bump;
          a2 = bField.getCompFrac(comp, cell, "pec");
          cell[comp] -= bump;

          aprods[j] = a1 * a2;
          std::cout << "    a_{" << comps[i / 2] << ((i % 2) ? '+' : '-')
                    << "} * a_{" << comps[j / 2] << ((j % 2) ? '+' : '-') 
                    << "} = " << aprods[j] << "\n";
          rowsum += aprods[j];

          if ((j >= i) and (aprods[j] / vfrac > large)) largeRatios++;
        }

        std::cout << "  ratio of area product sum to vfrac = " << rowsum / vfrac << "\n";
      }
    }
  }
  std::cout << "Number of a1*a2 / v greater than " << large << ": " << largeRatios << "\n";

}
#endif

template<size_t DIM, typename Scalar>
void MxMagWaveOp<DIM, Scalar>::checkEigensolution(
std::vector<std::complex<double> > const & eigVals,
MxMultiVector<Scalar> const & eigVecs,
std::vector<double> & residuals) const {


  MxMultiVector<Scalar> Aevecs(eigVecs.getMap(), eigVecs.getNumVecs());
  MxMultiVector<Scalar> Bevecs(eigVecs.getMap(), eigVecs.getNumVecs());

  // apply operators; here we'll apply the curl curl operator, not
  // the vector laplacian as we did in the solve
  //levelMats[0].find("vecLapl")->second->Apply(eigVecs, *Aevecs);
  mEMOps->getOp("curlCurl")->apply(eigVecs, Aevecs);
  mRhs->apply(eigVecs, Bevecs);

#if 0
  // allocations
  Epetra_MultiVector * Aevecs, * Bevecs, * bWork(0); // where Ax = \lambda Bx
  Aevecs = new Epetra_MultiVector(eigVecs.Map(), eigVals.size());
  if (mSim->hasPEC() and mNeedMu) {
    Bevecs = new Epetra_MultiVector(eigVecs.Map(), eigVals.size());
    bWork = new Epetra_MultiVector(eigVecs.Map(), eigVals.size());
  }
  else if (mSim->hasPEC() or mNeedMu)
    Bevecs = new Epetra_MultiVector(eigVecs.Map(), eigVals.size());
  else
    // no worries, won't be modified
    Bevecs = const_cast<Epetra_MultiVector *>(&eigVecs);

  Epetra_MultiVector * eWork = new Epetra_MultiVector(*eMap, eigVals.size());
  levelMats[0].find("curlB")->second->Apply(eigVecs, *eWork);
  if (mNeedEps)
    levelMats[0].find("invEps")->second->Apply(*eWork, *eWork);
  if (mSim->hasPEC())
    levelMats[0].find("dmL")->second->Apply(*eWork, *eWork);
  levelMats[0].find("curlE")->second->Apply(*eWork, *Aevecs);
  delete eWork;


  if (mSim->hasPEC() and mNeedMu) {
    levelMats[0].find("mu")->second->Apply(eigVecs, *bWork);
    levelMats[0].find("dmA")->second->Apply(*bWork, *Bevecs);
  }
  else if (mNeedMu)
    levelMats[0].find("mu")->second->Apply(eigVecs, *Bevecs);
  else if (mSim->hasPEC())
    levelMats[0].find("dmA")->second->Apply(eigVecs, *Bevecs);
#endif



  // get k^2 values from eigenvalues
  std::vector<std::complex<double> > k2; 
  for (size_t i = 0; i < eigVals.size(); ++i) {
    if (invert) {
      k2.push_back(1.0 / eigVals[i] + shift);
    }
    else {
      k2.push_back(eigVals[i] + shift);
    }
  }

  // form A x - k^2 B x = 0
  Scalar tmp;
  for (size_t i = 0; i < eigVals.size(); ++i) {
    // ignoring complex eigenvalues for now, they represent
    // asymmetries in the mag wave operator (because of dielectric
    // update). Residuals will be higher for these modes.
    //Aevecs(i)->Update(-k2[i].real(), *Bevecs(i), 1.0);
    // shallow copy ------|
    MxUtil::convertScalar(-k2[i], tmp);
    Aevecs.getVectorNonConst(i, false)->update(tmp, *Bevecs.getVector(i, false),
        ScalarTraits<Scalar>::one());
  }

  // get residuals
  Aevecs.norm2(residuals);

  // scale to get relative error
  for (size_t i = 0; i < k2.size(); ++i)
    residuals[i] /= abs(k2[i]);

#if 0
  delete Aevecs;
  if (mSim->hasPEC() and mNeedMu)
    delete Bevecs, bWork;
  else if (mSim->hasPEC() or mNeedMu)
    delete Bevecs;
#endif

}

template<size_t DIM, typename Scalar>
void MxMagWaveOp<DIM, Scalar>::checkDivergences(
MxMultiVector<Scalar> const & magVecs,
std::vector<double> & residuals) const {
  size_t nVecs = magVecs.getNumVecs();
  residuals.resize(nVecs);

#if 0
  Epetra_MultiVector b(magVecs);
  Epetra_MultiVector psi(*psiMap, nVecs, false);
  if (mNeedMu)
    levelMats[0].find("mu")->second->Apply(magVecs, b);
  if (mSim->hasPEC())
    levelMats[0].find("dmA")->second->Apply(b, b);

  levelMats[0].find("divB")->second->Apply(b, psi);
#endif

  //MxMultiVector<Scalar> b(magVecs), psi(mPsiMap, nVecs);
  MxMultiVector<Scalar> b(mBMap, nVecs), psi(mPsiMap, nVecs);
  mRhs->apply(magVecs, b);
  mEMOps->getOp("divB")->apply(b, psi);
  psi.norm2(residuals);
}


template<size_t DIM, typename Scalar>
void MxMagWaveOp<DIM, Scalar>::magToElec(
MxMultiVector<Scalar> const & mag,
MxMultiVector<Scalar> & elec) const {

#if 0
  levelMats[0].find("curlB")->second->Apply(mag, elec);
  if (mSim->hasDielectric())
    levelMats[0].find("invEps")->second->Apply(elec, elec);
#endif
  mEMOps->getOp("curlB")->apply(mag, elec);
  if (mSim->hasDielectric())
    mEMOps->getOp("invEps")->apply(elec, elec);
}

template<size_t DIM, typename Scalar>
void MxMagWaveOp<DIM, Scalar>::eigValsToFreqs(
std::vector<std::complex<double> > const & eigVals,
std::vector<std::complex<double> > & freqs) const {
  freqs.resize(eigVals.size());

  // get k^2 values from eigenvalues
  std::complex<double> val;
  for (size_t i = 0; i < eigVals.size(); ++i) {
    if (invert) {
      val = 1.0 / eigVals[i] + shift;
    }
    else {
      val = eigVals[i] + shift;
    }
    // output in Hz
    freqs[i] = sqrt(val) * MxUtil::lightspeed / 2.0 / MxUtil::pi;
  }
}



template class MxMagWaveOp<1, double>;
template class MxMagWaveOp<2, double>;
template class MxMagWaveOp<3, double>;
template class MxMagWaveOp<1, MxComplex>;
template class MxMagWaveOp<2, MxComplex>;
template class MxMagWaveOp<3, MxComplex>;
