
#include "MxWakeSolver.h"

#include <iostream>
#include <string>
#include <complex>

#include "MxGrid.h"
#include "MxGridField.hpp"
#include "MxMagWaveOp.h"
#include "MxIO.h"

#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"


template<size_t DIM>
MxWakeSolver<DIM>::MxWakeSolver(const MxEMSim<DIM> * aSim, Teuchos::ParameterList aPList) : 
myPID(aSim->getGrid().getComm().MyPID()), sim(aSim), pList(aPList) {
  MxMagWaveOp<DIM> * mwOp = new MxMagWaveOp<DIM>(sim);
  op = Teuchos::rcp(mwOp);

  double time = pList.get("wake solver : time", 0.0);
  double nyFreq = pList.get("wake solver : max frequency", 0.0);
  double numFreqs = pList.get("wake solver : num frequencies", 0.0);
  MxDimVector<double, 2>
  mSigma = pList.get("wake solver : bunch sigma", 0.0);
  mRadius = pList.get("wake solver : bunch radius", 0.0);
  mCenter = pList.get("wake solver : bunch center", MxDimVector<double, 2>(0));


  double dt = 0.5 / nyFreq;
  int steps = int(time / dt);


  if (time == 0.0) {
    exit(EXIT_FAILURE);
  }


  


  int blockSize = pList.get("eigensolver : block size", 1);
  if (mwOp->isComplex())
    blockSize *= 2;
  int numBlocks = pList.get("eigensolver : basis", 20);
  int maxRestarts = pList.get("eigensolver : max restarts", 20);
  double tol = pList.get("eigensolver : tol", 1.e-8);
  int stepSize = pList.get("eigensolver : step size", 1);
  std::string ortho = pList.get("eigensolver : ortho method", "SVQB");
  std::string targetSpectrum = pList.get("eigensolver : spectrum", "LM");
  int nev = pList.get("eigensolver : nev", 10);
  int verb = pList.get("eigensolver : output", 1);

  anasaziPList.set("Verbosity", verb);
  anasaziPList.set("Which", targetSpectrum);
  anasaziPList.set("Block Size", blockSize);
  anasaziPList.set("Num Blocks", numBlocks);
  anasaziPList.set("Maximum Restarts", maxRestarts);
  anasaziPList.set("Convergence Tolerance", tol);
  anasaziPList.set("Step Size", stepSize);
  anasaziPList.set("Orthogonalization", ortho);

  //initVec = Teuchos::rcp(new Epetra_MultiVector(sim->getField("bfield").getMap(), blockSize));

  initVec = Teuchos::rcp(new Epetra_MultiVector(op->OperatorDomainMap(), blockSize));
  initVec->Random();

  uniFields = Teuchos::rcp(new Epetra_MultiVector(sim->getField("bfield").uniformFields()));

  eigProb = Teuchos::rcp(new Anasazi::BasicEigenproblem<double, MV, OP>(op, initVec));
  eigProb->setNEV(nev);
  //if (!sim->hasDielectric())
    //eigProb->setHermitian(true);
  //eigProb->setAuxVecs(uniFields);

  bool boolret = eigProb->setProblem();
  if (boolret != true) {
    if (myPID == 0) {
      std::cout << "Anasazi::BasicEigenproblem::setProblem() \
                returned with error." << std::endl;
    }
  }

  eigMgr = Teuchos::rcp(new Anasazi::BlockKrylovSchurSolMgr<double, MV, OP>(eigProb, anasaziPList));

}


template<size_t DIM>
void MxWakeSolver<DIM>::solveOmega(double omega) {

  fillJ(omega, *mJ);
  mwOp->getOp("curlE")->Apply(*mJ, *mCurlJ);

  mSolver->SetRHS(mCurlJ.get());

  mwOp->Apply(*mCurlJ, *mB);

  // now get E from B and J

}

template<size_t DIM>
void MxWakeSolver<DIM>::fillJ(double omega, Epetra_MultiVector & j) {
  
  MxGridFieldIter<DIM> iter(&sim->getField("efield"));

  MxComplex val, i(0.0, 1.0);

  double coeff = omega * mSigma / (sqrt(2.0) * MxUtil::lightspeed);
  coeff = mCharge * MxUtil::lightspeed * exp(-coeff * coeff);

  int comp, row;
  MxDimVector<double, DIM> coord;
  double r;

  for (iter.begin(); !iter.atEnd(); iter.bump()) {
    comp = iter.getComp();
    coord = iter.getCoord();
    row = iter.getGlobCompIndx();

    // component in bunch?
    r = (MxDimVector<double, 2>(coord) - mCenter).norm();

    // only set Jz
    if (comp == 2 and r < mRadius) {
      val = coeff * exp(i * omega * coord[2] / MxUtil::lightspeed);

      j.ReplaceGlobalValue(2*row, 0, val.real());
      j.ReplaceGlobalValue(2*row + 1, 0, val.imag());
    }
    else {
      j.ReplaceGlobalValue(2*row, 0, 0);
      j.ReplaceGlobalValue(2*row + 1, 0, 0);
    }

  }

}

template<size_t DIM>
void MxWakeSolver<DIM>::solve() {
  Anasazi::ReturnType returnCode = eigMgr->solve();
  if (returnCode != Anasazi::Converged && myPID==0) {
    std::cout << "Anasazi::EigensolverMgr::solve() \
                  returned unconverged." << std::endl;
  }

  Anasazi::Eigensolution<double, MV> sol = eigProb->getSolution();
  eigVecs = sol.Evecs;
  eigVals.resize(sol.numVecs);
  for (int k = 0; k < sol.numVecs; k++)
    eigVals[k] = std::complex<double>(sol.Evals[k].realpart, sol.Evals[k].imagpart);


  //sim->getField("bfield").save(*eigVecs, "eigenvectors");

  // analyze
  MxMagWaveOp<DIM> * mwOp = dynamic_cast<MxMagWaveOp<DIM> *>(op.get());

  // check solution
  mwOp->checkEigensolution(eigVals, *eigVecs, eigErrs);
  mwOp->eigValsToFreqs(eigVals, eigFreqs);

  if (pList.get("has curl null", false)) {
    std::vector<double> residuals;
    mwOp->checkDivergences(*eigVecs, residuals);
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
    for (size_t i = 0; i < eigFreqs.size(); ++i) {
      std::cout << std::setw(3) << i << ":" 
        << std::setw(16) << eigFreqs[i].real()
        << std::setw(20) << eigFreqs[i].imag()
        << std::setw(20) << eigErrs[i] << "\n";
    }
    std::cout << "\n=========================================\n";
  }

  // ugh, separate real/imag parts...
  std::vector<double> freqsRe, freqsIm;
  for (size_t i = 0; i < sol.numVecs; ++i) {
    freqsRe.push_back(eigFreqs[i].real());
    freqsIm.push_back(eigFreqs[i].imag());
  }
  MxUtil::HDF5::saveArray(&freqsRe[0], freqsRe.size(), "mxEigenFreqsReal.h5");
  MxUtil::HDF5::saveArray(&freqsIm[0], freqsIm.size(), "mxEigenFreqsImag.h5");

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

  MxIO<DIM> io(&sim->getGrid().getComm());
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

}


template<size_t DIM>
Teuchos::RCP<Epetra_MultiVector> MxWakeSolver<DIM>::getElecEigVecs(
bool realPart) const {
  if (realPart && elecEvecsRe != Teuchos::null) {
    return elecEvecsRe;
  }
  else if (!realPart && elecEvecsIm != Teuchos::null) {
    return elecEvecsIm;
  }
  else {
    std::cout << "Tried to obtain electric field eigenvectors from\n"
              << "MxWakeSolver without success. Returning null.\n";
    return Teuchos::null;
  }
}

template<size_t DIM>
Teuchos::RCP<Epetra_MultiVector> MxWakeSolver<DIM>::getMagEigVecs(
bool realPart) const {
  if (realPart && magEvecsRe != Teuchos::null) {
    return magEvecsRe;
  }
  else if (!realPart && magEvecsIm != Teuchos::null) {
    return magEvecsIm;
  }
  else {
    std::cout << "Tried to obtain magnetic field eigenvectors from\n"
              << "MxWakeSolver without success. Returning null.\n";
    return Teuchos::null;
  }
}

template<size_t DIM>
void MxWakeSolver<DIM>::saveFieldValsAtPoints(std::vector<MxDimVector<double, DIM> > const & points) const {
  MxIO<DIM> io(&sim->getGrid().getComm());
  io.save(points, *magEvecsRe, sim->getField("bfield"), "eigenvectors_at_points");
  io.save(points, *elecEvecsRe, sim->getField("efield"), "eigenvectors_at_points");
}

template<size_t DIM>
void MxWakeSolver<DIM>::saveFieldValsAtPoints(MxPointCloud<DIM> const & pointCloud) const {
  MxIO<DIM> io(&sim->getGrid().getComm());
  io.save(pointCloud, *magEvecsRe, sim->getField("bfield"), "eigenvectors_at_points");
  io.save(pointCloud, *elecEvecsRe, sim->getField("efield"), "eigenvectors_at_points");
}


template class MxWakeSolver<1>;
template class MxWakeSolver<2>;
template class MxWakeSolver<3>;
