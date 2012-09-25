#include "MxEMOps.h"

// Mx operators! yey!
#include "MxYeeDeyMittraCurlE.h"
#include "MxYeeDeyMittraCurlB.h"
#include "MxYeeDeyMittraDivB.h"
#include "MxYeeDeyMittraGradPsi.h"
#include "MxYeeDeyMittraFracs.h"
#include "MxYeeFitInvEps.h"
#include "MxYeeFitMu.h"
#include "MxEikxOp.h"

template<size_t DIM, typename Scalar>
MxEMOps<DIM, Scalar>::MxEMOps(RCP<MxEMSim<DIM> > sim,
bool initSingleOps) : 
pid(sim->getGrid().getComm()->myPID()), mSim(sim) {
  if (initSingleOps)
    setSingleOps();
}


template<size_t DIM, typename Scalar>
void MxEMOps<DIM, Scalar>::printRCPs() {
  typename std::map<std::string, RCP<MxCrsMatrix<Scalar> > >::const_iterator iter;
  for (iter = mOps.begin(); iter != mOps.end(); iter++) {
    if (pid == 0) std::cout << iter->first << "\n"
        << "  strong: " << iter->second.strong_count() << "\n"
        << "  weak: " << iter->second.weak_count() << "\n"
        << "  has own: " << iter->second.has_ownership() << "\n"
        << "  address: " << iter->second.get() << "\n";
  }
}

template<size_t DIM, typename Scalar>
void MxEMOps<DIM, Scalar>::setSingleOps(RCP<MxEMSim<DIM> > sim) {
  Teuchos::ParameterList pList = sim->getParameters();
  MxPolType pol = pList.get("polarization", TM);
  bool hasCurlNull = ((DIM == 3) or (DIM == 2 and pol == TM));

  const MxGrid<DIM> & grid = sim->getGrid();

  if (pid == 0) std::cout << "Constructing curl of E...\n";
  mOps.insert(std::make_pair("curlE", rcp(
    new MxYeeDeyMittraCurlE<DIM, Scalar>(sim))));
  if (pid == 0) std::cout << "Constructing curl of B...\n";
  mOps.insert(std::make_pair("curlB", rcp(
    new MxYeeDeyMittraCurlB<DIM, Scalar>(sim))));


  if (sim->hasPEC()) {
    if (pid == 0) std::cout << "Constructing Dey-Mittra area fracs...\n";
    mOps.insert(std::make_pair("dmA", rcp(
      new MxYeeDeyMittraFracs<DIM, Scalar>(
        &sim->getField("bfield"), &grid, false, 0.e-12))));

    if (pid == 0) std::cout << "Constructing Dey-Mittra length fracs...\n";
    mOps.insert(std::make_pair("dmL", rcp(
      new MxYeeDeyMittraFracs<DIM, Scalar>(
        &sim->getField("efield"), &grid, false, 0.e-6))));
  }

  if (sim->hasDielectric() or sim->hasPML()) {
    //symmetrizeInvEps();
#if 1
    if (pid == 0)  std::cout << "Constructing FIT inverse epsilon operator...\n";

    // asymmetric accurate algorithm
    MxYeeFitInvEps<DIM, Scalar> * invEps =
        new MxYeeFitInvEps<DIM, Scalar>(sim);

    //EpetraExt::RowMatrixToMatrixMarketFile("MxInvEps.dat", *invEps);

    mOps.insert(std::make_pair("invEps", rcp(invEps)));
    if (hasCurlNull)
      mOps.insert(std::make_pair("invEpsVolAve",
        invEps->cellAveInvEpsOperator()));
#endif

#if 0
    if (pid == 0)  std::cout << "Constructing MIT inverse epsilon operator...\n";

    // MIT symmetric inaccurate algorithm
    MxYeeMitInvEps<DIM> * invEps = new MxYeeMitInvEps<DIM>(&sim);
    //EpetraExt::RowMatrixToMatrixMarketFile("MxInvEps.dat", *invEps);
    mOps.insert(std::make_pair("invEps", Teuchos::rcp(invEps)));
    if (DIM == 3) {
      Epetra_CrsMatrix * invEpsVolAve = new Epetra_CrsMatrix(invEps->cellAveInvEpsOperator());
      mOps.insert(std::make_pair("invEpsVolAve", Teuchos::rcp(invEpsVolAve)));
    }
#endif

#if 0
    // symmetrized by 0.5 * (A + A^T)
    MxYeeFitInvEps<DIM, Scalar> * invEps = new MxYeeFitInvEps<DIM, Scalar>(sim);
    MxCrsMatrix<Scalar> * invEpsSymm = new MxCrsMatrix<Scalar>(invEps->getRangeMap());

    MxCrsMatrix<Scalar>::matrixMatrixAdd(*invEps, true, 0.5, *invEps, false, 0.5, *invEpsSymm, true);

    mOps.insert(std::make_pair("invEps", rcp(invEpsSymm)));

    if (hasCurlNull)
      mOps.insert(std::make_pair("invEpsVolAve",
        invEps->cellAveInvEpsOperator()));

    delete invEps;
#endif
  }

  if (hasCurlNull) {
    if (pid == 0) std::cout << "Constructing divergence of B...\n";
    mOps.insert(std::make_pair("divB", rcp(
      new MxYeeDeyMittraDivB<DIM, Scalar>(sim))));
    if (pid == 0) std::cout << "Constructing gradient of Psi...\n";
    mOps.insert(std::make_pair("gradPsi", rcp(
      new MxYeeDeyMittraGradPsi<DIM, Scalar>(sim))));

    if (sim->hasPEC()) {
      if (pid == 0) std::cout << "Constructing inverse Dey-Mittra volume fracs...\n";
      //conditionDeyMittraVolumes(*dmVInv, *mOps.find("dmA")->second, sim);
      mOps.insert(std::make_pair("dmVInv", rcp(
        new MxYeeDeyMittraFracs<DIM, Scalar>(
          &sim->getField("psifield"), &grid, true, 0.e-6, false, 0.3))));
    }
  }

  if (sim->hasMu() or sim->hasPML()) {
    if (pid == 0)  std::cout << "Constructing FIT inverse mu operator...\n";

    // asymmetric accurate algorithm
    MxYeeFitMu<DIM, Scalar> * mu =
      new MxYeeFitMu<DIM, Scalar>(sim, false); // not inverted

    //EpetraExt::RowMatrixToMatrixMarketFile("MxInvMu.dat", *mu);

    mOps.insert(std::make_pair("mu", rcp(mu)));
    if (hasCurlNull)
      mOps.insert(std::make_pair("invMuVolAve",
        mu->cellAveInvMuOperator()));
  }


#if 0
  if (mIsComplex) {
    MxDimVector<double, DIM> blochK(pList.get("phase shifts", MxDimVector<double, DIM>(0)));
    std::cout << "blochK : " << blochK << "\n";
    blochK /= grid.getSize();

    double phShift = 0.0;
    MxEikxOp<DIM> eikxMagRe(blochK, 1.0, phShift*MxUtil::pi, &aSim.getField("bfield"), &grid, false);
    MxEikxOp<DIM> eikxMagIm(blochK, 1.0, phShift*MxUtil::pi, &aSim.getField("bfield"), &grid, true);
    insertComplexOp("eikxMag", eikxMagRe, eikxMagIm, ops);

    if (mHasCurlNull) {
      MxEikxOp<DIM> eikxPsiRe(blochK, 1.0, phShift*MxUtil::pi,
          &aSim.getField("psifield"), &grid, false);
      MxEikxOp<DIM> eikxPsiIm(blochK, 1.0, phShift*MxUtil::pi,
          &aSim.getField("psifield"), &grid, true);
      insertComplexOp("eikxPsi", eikxPsiRe, eikxPsiIm, ops);
    }
  }
#endif

}

template<size_t DIM, typename Scalar>
RCP<MxCrsMatrix<Scalar> > MxEMOps<DIM, Scalar>::getOp(std::string name) const {

  typename std::map<std::string, RCP<MxCrsMatrix<Scalar> > >::const_iterator it;
  it = mOps.find(name);

  if (it == mOps.end()) {
    std::cout << "MxEMOps::getOp: operator '" << name
      << "' not found.\n";
    exit(EXIT_FAILURE);
  }

  return it->second;
}

template<size_t DIM, typename Scalar>
void MxEMOps<DIM, Scalar>::gatherOps(std::vector<std::string> const & names,
std::vector<RCP<MxCrsMatrix<Scalar> > > & mats) const {
  mats.clear();

  typename std::map<std::string, RCP<MxCrsMatrix<Scalar> > >::const_iterator it;
  for (size_t i = 0; i < names.size(); ++i)
    // checks existence
    mats.push_back(getOp(names[i]));
}

template<size_t DIM, typename Scalar>
RCP<MxCrsMatrix<Scalar> > MxEMOps<DIM, Scalar>::multiplyOps(
std::vector<std::string> names,
std::vector<bool> transposes, bool fillComplete, bool store,
std::string name) {
  
  std::vector<RCP<MxCrsMatrix<Scalar> > > mats;
  gatherOps(names, mats);

  RCP<MxMap> rowMap, colMap;
  if (transposes.back())
    rowMap = mats.back()->getDomainMap();
  else
    rowMap = mats.back()->getRangeMap();
    
  if (transposes[0])
    colMap = mats[0]->getRangeMap();
  else
    colMap = mats[0]->getDomainMap();
  
  //RCP<MxCrsMatrix<Scalar> > res = rcp(new MxCrsMatrix<Scalar>(rowMap, colMap));
  RCP<MxCrsMatrix<Scalar> > res = rcp(new MxCrsMatrix<Scalar>(rowMap));

  MxUtil::Trilinos::massiveCrsMultiply(mats, transposes, res, fillComplete);

  if (store)
    mOps.insert(std::make_pair(name, res));
  
  return res;
}


template<size_t DIM, typename Scalar>
RCP<MxCrsMatrix<Scalar> > MxEMOps<DIM, Scalar>::addOps(
std::vector<std::string> names,
std::vector<Scalar> coeffs, bool store,
std::string name) {

  std::vector<RCP<MxCrsMatrix<Scalar> > > mats;
  gatherOps(names, mats);

  RCP<MxMap> rowMap = mats[0]->getRangeMap();
  RCP<MxMap> colMap = mats[0]->getDomainMap();

  RCP<MxCrsMatrix<Scalar> > m1, m2, res;
  Scalar s1, s2;

  for (size_t i = 1; i < mats.size(); ++i) {
    if (i == 1) {
      m1 = mats[0];
      m2 = mats[1];
      s1 = coeffs[0];
      s2 = coeffs[1];
    }
    else {
      m1 = res;
      m2 = mats[i];
      s1 = ScalarTraits<Scalar>::one();
      s2 = coeffs[i];
    }
    //res = rcp(new MxCrsMatrix<Scalar>(rowMap, colMap));
    res = rcp(new MxCrsMatrix<Scalar>(rowMap));

    MxCrsMatrix<Scalar>::matrixMatrixAdd(
      *m1, false, s1, *m2, false, s2, *res, true);
  }

  if (store)
    mOps.insert(std::make_pair(name, res));
  
  return res;

}


template<size_t DIM, typename Scalar>
void MxEMOps<DIM, Scalar>::saveOps() const {
  if (pid == 0) std::cout << "Saving operators in MxEMOps...\n";
  typename std::map<std::string, RCP<MxCrsMatrix<Scalar> > >::const_iterator iter;
  for (iter = mOps.begin(); iter != mOps.end(); iter++)
    iter->second->save(iter->first);
}


template class MxEMOps<1, double>;
template class MxEMOps<2, double>;
template class MxEMOps<3, double>;
template class MxEMOps<1, MxComplex>;
template class MxEMOps<2, MxComplex>;
template class MxEMOps<3, MxComplex>;
