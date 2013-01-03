#include "MxSolver.h"

#include <iostream>
#include <string>
#include <complex>

#include "MxGrid.h"
#include "MxGridField.hpp"
#include "MxMagWaveOp.h"
#include "MxIO.h"
#include "MxAnasaziMV.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"

#include "Epetra_Comm.h"
#include "Epetra_Time.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"


template<size_t DIM, typename Scalar>
MxSolver<DIM, Scalar>::MxSolver(RCP<MxEMSim<DIM> > aSim, Teuchos::ParameterList aPList) : 
myPID(aSim->getGrid().getComm()->myPID()), sim(aSim), pList(aPList) {
  MxMagWaveOp<DIM, Scalar> * mwOp = new MxMagWaveOp<DIM, Scalar>(sim);
  op = rcp(mwOp);

  std::string type = pList.get("eigensolver : type", "krylov-schur");
  int blockSize = pList.get("eigensolver : block size", 1);
  if (pList.get("is complex", false))
    blockSize *= 2;
  int numBlocks = pList.get("eigensolver : basis", 20);
  int maxRestarts = pList.get("eigensolver : max restarts", 20);
  double tol = pList.get("eigensolver : tol", 1.e-8);
  int stepSize = pList.get("eigensolver : step size", 1);
  std::string ortho = pList.get("eigensolver : ortho method", "SVQB");
  std::string targetSpectrum = pList.get("eigensolver : spectrum", "LM");
  int nev = pList.get("eigensolver : nev", 10);
  int verb = pList.get("eigensolver : output", 1);

  //anasaziPList.set("Verbosity", verb);
  anasaziPList.set("Verbosity", 10);
  anasaziPList.set("Which", targetSpectrum);
  anasaziPList.set("Block Size", blockSize);
  anasaziPList.set("Num Blocks", numBlocks);
  anasaziPList.set("Maximum Restarts", maxRestarts);
  anasaziPList.set("Convergence Tolerance", tol);
  anasaziPList.set("Step Size", stepSize);
  anasaziPList.set("Orthogonalization", ortho);

  std::cout << anasaziPList;


  // epetra interface...
#if 0
  mInitVec = rcp(new Epetra_MultiVector(op->OperatorDomainMap(), blockSize));
  mInitVec->Random();
  mEigProb = rcp(new Anasazi::BasicEigenproblem<double, MV, OP>(op, mInitVec));
#endif

// my Anasazi interface
#if 1
  RCP<MxAnasaziMV<Scalar> > initVec =
    rcp(new MxAnasaziMV<Scalar>(sim->getField("bfield").getMap(), blockSize));
  initVec->random();
  MxGridField<DIM>::zeroUnusedComponents(*initVec, sim->getField("bfield"));
  mInitVec = initVec;

  mEigProb = rcp(new Anasazi::BasicEigenproblem<Scalar, MV, OP>(op, mInitVec));
#endif

  mEigProb->setNEV(nev);
  if (!sim->hasDielectric() and !sim->hasMu())
    mEigProb->setHermitian(true);

  //mEigProb->setAuxVecs(uniFields);

  bool boolret = mEigProb->setProblem();
  if (boolret != true) {
    if (myPID == 0) {
      std::cout << "Anasazi::BasicEigenproblem::setProblem() \
                returned with error." << std::endl;
    }
  }

  if (type == "krylov-schur") {
    mEigMgr = rcp(new Anasazi::BlockKrylovSchurSolMgr<Scalar, MV, OP>(mEigProb, anasaziPList));
  } else if (type == "davidson") {
    anasaziPList.set("Use Locking", true);
    mEigMgr = rcp(new Anasazi::BlockDavidsonSolMgr<Scalar, MV, OP>(mEigProb, anasaziPList));
  } else {
    std::cout << "Unsupported eigensolver type: " << type << "\n";
    exit(EXIT_FAILURE);
  //mEigMgr = rcp(new Anasazi::BlockKrylovSchurSolMgr<double, MV, OP>(mEigProb, anasaziPList));
  }

}


template<size_t DIM, typename Scalar>
void MxSolver<DIM, Scalar>::solve() {
  
  Epetra_Time timer(*sim->getGrid().getComm()->getEpetraComm());
  Anasazi::ReturnType returnCode = mEigMgr->solve();
  if (returnCode != Anasazi::Converged && myPID==0) {
    std::cout << "Anasazi::EigensolverMgr::solve() \
                  returned unconverged." << std::endl;
  }
  if (myPID == 0) {
    std::cout << "\n\n"
              << "---------------------------\n"
              << " Eigensolve time: " << timer.ElapsedTime() << " seconds\n"
              << "---------------------------\n\n";
  }


#if 1
  Anasazi::Eigensolution<Scalar, MV> sol = mEigProb->getSolution();
  mBEigVecs = Teuchos::rcp_dynamic_cast<MxAnasaziMV<Scalar> >(sol.Evecs, true); // throw on fail
#else
  Anasazi::Eigensolution<double, MV> sol = mEigProb->getSolution();
  mBEigVecs = rcp(new MxMultiVector<Scalar>(sim->getField("bfield").getMap(), sol.numVecs));
  mBEigVecs->setRawMV(sol.Evecs);
#endif

  mEigVals.resize(sol.numVecs);
  for (int k = 0; k < sol.numVecs; k++)
    mEigVals[k] = MxComplex(sol.Evals[k].realpart, sol.Evals[k].imagpart);

#if 0
  // Print norms, max values
  std::vector<double> norms(sol.numVecs);
  mBEigVecs->norm2(norms);
  Scalar val;
  MxDimVector<int, DIM> cell;
  int comp;
  std::cout << "\n\n========== Max Values ==========\n";
  for (int i = 0; i < sol.numVecs; ++i) {
    sim->getField("bfield").maxValueLocation(*mBEigVecs, i, cell, comp, val);
  }
  std::cout << "===========================\n\n";
#endif

  // analyze
  RCP<MxMagWaveOp<DIM, Scalar> > mwOp = Teuchos::rcp_dynamic_cast<MxMagWaveOp<DIM, Scalar> >(op, true);

  // check solution
  mwOp->checkEigensolution(mEigVals, *mBEigVecs, eigErrs);
  mwOp->eigValsToFreqs(mEigVals, mEigFreqs);

  if (pList.get("has curl null", false)) {
    std::vector<double> residuals;
    mwOp->checkDivergences(*mBEigVecs, residuals);
    MxUtil::printStdVector(residuals);
  }

  // print and save solution info
  if (myPID == 0) { 
    std::cout << "\n========= Eigenfrequencies (Hz) =========\n\n";
    std::cout.setf(std::ios_base::left, std::ios_base::adjustfield);
    std::cout << std::setw(20) << "Real"
              << std::setw(20) << "Imag"
              << std::setw(20) << "Relative error" << "\n\n";

    std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);
    std::cout.setf(std::ios::scientific);
    std::cout.precision(8);
    for (size_t i = 0; i < mEigFreqs.size(); ++i) {
      std::cout << std::setw(3) << i << ":" 
        << std::setw(16) << mEigFreqs[i].real()
        << std::setw(20) << mEigFreqs[i].imag()
        << std::setw(20) << eigErrs[i] << "\n";
    }
    std::cout << "\n=========================================\n";
  }

  // ugh, separate real/imag parts...
  std::vector<double> freqsRe, freqsIm;
  for (size_t i = 0; i < sol.numVecs; ++i) {
    freqsRe.push_back(mEigFreqs[i].real());
    freqsIm.push_back(mEigFreqs[i].imag());
  }
  //MxUtil::HDF5::saveArray(&freqsRe[0], freqsRe.size(), "mxEigenFreqsReal.h5");
  //MxUtil::HDF5::saveArray(&freqsIm[0], freqsIm.size(), "mxEigenFreqsImag.h5");

  // get/save magnetic and electric fields
  mEEigVecs = Teuchos::rcp(new MxMultiVector<Scalar>(
      sim->getField("efield").getMap(), sol.numVecs));
  mwOp->magToElec(*mBEigVecs, *mEEigVecs);

#if 0
  // get/save magnetic and electric fields
  if (mwOp->isComplex()) {
    // allocate
    elecEvecs = Teuchos::rcp(new Epetra_MultiVector(mwOp->getEMap(), sol.numVecs));
    magEvecsRe = Teuchos::rcp(new Epetra_MultiVector(
        sim->getField("bfield").getMap(), sol.numVecs));
    magEvecsIm = Teuchos::rcp(new Epetra_MultiVector(
        sim->getField("bfield").getMap(), sol.numVecs));
    elecEvecsRe = Teuchos::rcp(new Epetra_MultiVector(
        sim->getField("efield").getMap(), sol.numVecs));
    elecEvecsIm = Teuchos::rcp(new Epetra_MultiVector(
        sim->getField("efield").getMap(), sol.numVecs));

    // get re/im parts of eigensolution
    mwOp->magToElec(*eigVecs, *elecEvecs);
    mwOp->splitComplexMV("mag", *eigVecs, *magEvecsRe, *magEvecsIm);
    mwOp->splitComplexMV("elec", *elecEvecs, *elecEvecsRe, *elecEvecsIm);
    //mwOp->splitComplexMV("mag", *eigVecs, *magEvecsRe, *magEvecsIm);
    //mwOp->magToElec(*magEvecsRe, *magEvecsIm, *elecEvecsRe, *elecEvecsIm);

  }
  else {
    // allocate (only need electric fields)
    elecEvecsRe = Teuchos::rcp(new Epetra_MultiVector(
        sim->getField("efield").getMap(), sol.numVecs));

    magEvecsRe = eigVecs;
    mwOp->magToElec(*magEvecsRe, *elecEvecsRe);
  }
#endif

  MxIO<DIM, Scalar> io(sim->getGrid().getComm());
  io.save(*mBEigVecs, sim->getField("bfield"), "eigVecs");
  io.save(*mEEigVecs, sim->getField("efield"), "eigVecs");

#if 0
  if (mwOp->isComplex()) {
    io.save(*magEvecsRe, sim->getField("bfield"), "eigVecsRe");
    io.save(*elecEvecsRe, sim->getField("efield"), "eigVecsRe");
    io.save(*magEvecsIm, sim->getField("bfield"), "eigVecsIm");
    io.save(*elecEvecsIm, sim->getField("efield"), "eigVecsIm");
  }
  else {
    io.save(*magEvecsRe, sim->getField("bfield"), "eigVecs");
    io.save(*elecEvecsRe, sim->getField("efield"), "eigVecs");
  }
#endif

}

#if 0

template<size_t DIM>
Teuchos::RCP<Epetra_MultiVector> MxSolver<DIM>::getElecEigVecs(
bool realPart) const {
  if (realPart && elecEvecsRe != Teuchos::null) {
    return elecEvecsRe;
  }
  else if (!realPart && elecEvecsIm != Teuchos::null) {
    return elecEvecsIm;
  }
  else {
    std::cout << "Tried to obtain electric field eigenvectors from\n"
              << "MxSolver without success. Returning null.\n";
    return Teuchos::null;
  }
}

template<size_t DIM>
Teuchos::RCP<Epetra_MultiVector> MxSolver<DIM>::getMagEigVecs(
bool realPart) const {
  if (realPart && magEvecsRe != Teuchos::null) {
    return magEvecsRe;
  }
  else if (!realPart && magEvecsIm != Teuchos::null) {
    return magEvecsIm;
  }
  else {
    std::cout << "Tried to obtain magnetic field eigenvectors from\n"
              << "MxSolver without success. Returning null.\n";
    return Teuchos::null;
  }
}

template<size_t DIM>
void MxSolver<DIM>::saveFieldValsAtPoints(std::vector<MxDimVector<double, DIM> > const & points) const {
  MxIO<DIM> io(&sim->getGrid().getComm());
  io.save(points, *magEvecsRe, sim->getField("bfield"), "eigenvectors_at_points");
  io.save(points, *elecEvecsRe, sim->getField("efield"), "eigenvectors_at_points");
}

template<size_t DIM>
void MxSolver<DIM>::saveFieldValsAtPoints(MxPointCloud<DIM> const & pointCloud) const {
  MxIO<DIM> io(&sim->getGrid().getComm());
  io.save(pointCloud, *magEvecsRe, sim->getField("bfield"), "eigenvectors_at_points");
  io.save(pointCloud, *elecEvecsRe, sim->getField("efield"), "eigenvectors_at_points");
}

#endif

template class MxSolver<1, double>;
template class MxSolver<2, double>;
template class MxSolver<3, double>;
template class MxSolver<1, MxComplex>;
template class MxSolver<2, MxComplex>;
template class MxSolver<3, MxComplex>;
