
#include "MxYeeFitMu.h"

#include <limits>

#include "MxDimVector.hpp"
#include "MxGridFieldIter.hpp"



template<size_t DIM, typename Scalar>
MxYeeFitMu<DIM, Scalar>::MxYeeFitMu(RCP<MxEMSim<DIM> > theSim,
bool invert) :
mInvert(invert),
//MxCrsMatrix<Scalar>(theSim->getField("hfield").getMap(),
//  theSim->getField("hfield").getMap()),
MxCrsMatrix<Scalar>(theSim->getField("hfield").getMap()),
sim(theSim), hfield(&theSim->getField("hfield")),
bfield(&theSim->getField("bfield")),
grid(&theSim->getGrid()) {

  if (hfield->getNumComps() == 1)
    setMatrix2dTE();
  else
    setMatrix();
}

template<size_t DIM, typename Scalar>
typename MxYeeFitMu<DIM, Scalar>::Stencil
MxYeeFitMu<DIM, Scalar>::getStencil(size_t comp0,
MxDimVector<int, DIM> cell) const {
  size_t comp1 = (comp0 + 1) % hfield->getNumComps();
  size_t comp2 = (comp1 + 1) % hfield->getNumComps();

  size_t stenSize = 1 + 4*(DIM - 1);

  typename MxYeeFitMu<DIM, Scalar>::Stencil stencil;
  stencil.values = std::vector<Scalar>(stenSize,
    ScalarTraits<Scalar>::zero());

  stencil.comps = std::vector<size_t>(stenSize);
  stencil.comps[0] = comp0;
  for (size_t i = 1; i < 5; ++i)
    stencil.comps[i] = comp1;
  if (DIM == 3)
    for (size_t i = 5; i < 9; ++i)
      stencil.comps[i] = comp2;

  stencil.cells = std::vector<MxDimVector<int, DIM> >(stenSize, cell);
  stencil.cells[1][comp0]--;
  stencil.cells[2][comp0]--; stencil.cells[2][comp1]++;
  stencil.cells[4][comp1]++;
  if (DIM == 3) {
    stencil.cells[5][comp0]--;
    stencil.cells[6][comp0]--; stencil.cells[6][comp2]++;
    stencil.cells[8][comp2]++;
  }

  stencil.columns = std::vector<MxIndex>(stenSize);
  for (size_t i = 0; i < stenSize; ++i)
    stencil.columns[i] =
      bfield->globCompIndx(stencil.comps[i], stencil.cells[i]);

  stencil.bcfactors = std::vector<MxComplex>(stenSize);
  for (size_t i = 0; i < stenSize; ++i)
    stencil.bcfactors[i] =
      bfield->getCompFactor(stencil.comps[i], stencil.cells[i]);

  return stencil;
}

#if 0
typedef MxYeeFitMu<2>::Stencil Stencil2D;
typedef MxYeeFitMu<3>::Stencil Stencil3D;

template<>
MxYeeFitMu<1>::Stencil MxYeeFitMu<1>::getStencil(int comp0,
MxDimVector<int, 1> cell) const {
  return MxYeeFitMu<1>::Stencil();
}

template<>
Stencil2D MxYeeFitMu<2>::getStencil(int comp0,
MxDimVector<int, 2> cell) const {
  int comp1 = (comp0 + 1) % 2;

  Stencil2D stencil;
  stencil.cvalues = std::vector<MxComplex>(5, 0.0);
  stencil.rvalues = std::vector<double>(5, 0.0);

  stencil.comps = std::vector<int>(5);
  stencil.comps[0] = comp0;
  for (size_t i = 1; i < 5; ++i)
    stencil.comps[i] = comp1;
  
  stencil.cells = std::vector<MxDimVector<int, 2> >(5, cell);
  stencil.cells[1][comp0]--;
  stencil.cells[2][comp0]--; stencil.cells[2][comp1]++;
  stencil.cells[4][comp1]++;

  stencil.columns = std::vector<int>(5);
  for (size_t i = 0; i < 5; ++i)
    stencil.columns[i] = bfield->globCompIndx(stencil.comps[i], stencil.cells[i]);

  stencil.bcfactors = std::vector<MxComplex>(5);
  for (size_t i = 0; i < 5; ++i)
    stencil.bcfactors[i] = bfield->getCompFactor(stencil.comps[i], stencil.cells[i]);

  return stencil;
}

template<>
Stencil3D MxYeeFitMu<3>::getStencil(int comp0,
MxDimVector<int, 3> cell) const {
  int comp1 = (comp0 + 1) % 3;
  int comp2 = (comp1 + 1) % 3;

  Stencil3D stencil;
  stencil.cvalues = std::vector<MxComplex>(9, 0.0);
  stencil.rvalues = std::vector<double>(9, 0.0);

  stencil.comps = std::vector<int>(9);
  stencil.comps[0] = comp0;
  for (size_t i = 1; i < 5; ++i)
    stencil.comps[i] = comp1;
  for (size_t i = 5; i < 9; ++i)
    stencil.comps[i] = comp2;
  
  stencil.cells = std::vector<MxDimVector<int, 3> >(9, cell);
  stencil.cells[1][comp0]--;
  stencil.cells[2][comp0]--; stencil.cells[2][comp1]++;
  stencil.cells[4][comp1]++;
  stencil.cells[5][comp0]--;
  stencil.cells[6][comp0]--; stencil.cells[6][comp2]++;
  stencil.cells[8][comp2]++;

  stencil.columns = std::vector<int>(9);
  for (size_t i = 0; i < 9; ++i)
    stencil.columns[i] = bfield->globCompIndx(stencil.comps[i], stencil.cells[i]);

  stencil.bcfactors = std::vector<MxComplex>(9);
  for (size_t i = 0; i < 9; ++i)
    stencil.bcfactors[i] = bfield->getCompFactor(stencil.comps[i], stencil.cells[i]);

  return stencil;
}
#endif


template<size_t DIM, typename Scalar>
void MxYeeFitMu<DIM, Scalar>::getTuples(
typename MxYeeFitMu<DIM, Scalar>::Stencil const & stencil,
std::vector<typename MxYeeFitMu<DIM, Scalar>::Tuple> & tuples) const {
  tuples.clear();

  const size_t numTuples = 4*(DIM - 1); // 4 for 2d, 8 for 3d

  // stencil indices
  size_t stencilInds[numTuples][DIM];
  if (DIM == 2) {
    stencilInds[0][0] = 0; stencilInds[0][1] = 1;
    stencilInds[1][0] = 0; stencilInds[1][1] = 2;
    stencilInds[2][0] = 0; stencilInds[2][1] = 3;
    stencilInds[3][0] = 0; stencilInds[3][1] = 4;
  }
  else if (DIM == 3) {
    stencilInds[0][0] = 0; stencilInds[0][1] = 1; stencilInds[0][2] = 5;
    stencilInds[1][0] = 0; stencilInds[1][1] = 1; stencilInds[1][2] = 6;
    stencilInds[2][0] = 0; stencilInds[2][1] = 2; stencilInds[2][2] = 5;
    stencilInds[3][0] = 0; stencilInds[3][1] = 2; stencilInds[3][2] = 6;
    stencilInds[4][0] = 0; stencilInds[4][1] = 3; stencilInds[4][2] = 7;
    stencilInds[5][0] = 0; stencilInds[5][1] = 3; stencilInds[5][2] = 8;
    stencilInds[6][0] = 0; stencilInds[6][1] = 4; stencilInds[6][2] = 7;
    stencilInds[7][0] = 0; stencilInds[7][1] = 4; stencilInds[7][2] = 8;
  }
  
  typename MxYeeFitMu<DIM, Scalar>::Tuple tuple;

  size_t stenInd;
  double pecFrac;
  bool use;
  for (size_t i = 0; i < numTuples; ++i) {

    use = true;

    for (size_t j = 0; j < DIM; ++j) {
      stenInd = stencilInds[i][j];
      tuple.stencilInds[j] = stenInd;
      tuple.comps[j] = stencil.comps[stenInd];
      tuple.cells[j] = stencil.cells[stenInd];
      tuple.columns[j] = stencil.columns[stenInd];

      if (sim->hasPEC()) {
        pecFrac = bfield->getCompFrac(tuple.comps[j], tuple.cells[j], "pec");
        if (pecFrac == 0)
          use = false;
      }
    }

    if (use)
      tuples.push_back(tuple);
  }
}


#if 0
typedef MxYeeFitMu<2>::Tuple Doublet;
typedef MxYeeFitMu<3>::Tuple Triplet;

template<>
void MxYeeFitMu<1>::getTuples(MxYeeFitMu<1>::Stencil const & stencil,
std::vector<MxYeeFitMu<1>::Tuple> & singlets) const {}

template<>
void MxYeeFitMu<2>::getTuples(Stencil2D const & stencil,
std::vector<Doublet> & doublets) const {
  doublets.clear();

  // doublet stencil indices
  int doubletStencilInds[4][2] = {{0, 1}, {0, 2}, {0, 3}, {0, 4}};
  
  Doublet doublet;

  int stenInd;
  double pecFrac;
  bool use;
  for (size_t i = 0; i < 4; ++i) {

    use = true;

    for (size_t j = 0; j < 2; ++j) {
      stenInd = doubletStencilInds[i][j];
      doublet.stencilInds[j] = stenInd;
      doublet.comps[j] = stencil.comps[stenInd];
      doublet.cells[j] = stencil.cells[stenInd];
      doublet.columns[j] = stencil.columns[stenInd];

      if (sim->hasPEC()) {
        pecFrac = bfield->getCompFrac(doublet.comps[j], doublet.cells[j], "pec");
        if (pecFrac == 0)
          use = false;
      }
    }

    if (use)
      doublets.push_back(doublet);
  }
}

template<>
void MxYeeFitMu<3>::getTuples(Stencil3D const & stencil,
std::vector<Triplet> & triplets) const {
  triplets.clear();

  // triplet stencil indices
  int tripStenInds[8][3] = {{0, 1, 5},
                            {0, 1, 6},
                            {0, 2, 5},
                            {0, 2, 6},
                            {0, 3, 7},
                            {0, 3, 8},
                            {0, 4, 7},
                            {0, 4, 8}};
  
  Triplet triplet;

  int stenInd;
  double pecFrac;
  bool use;
  for (size_t i = 0; i < 8; ++i) {

    use = true;

    for (size_t j = 0; j < 3; ++j) {
      stenInd = tripStenInds[i][j];
      triplet.stencilInds[j] = stenInd;
      triplet.comps[j] = stencil.comps[stenInd];
      triplet.cells[j] = stencil.cells[stenInd];
      triplet.columns[j] = stencil.columns[stenInd];

      if (sim->hasPEC()) {
        pecFrac = hfield->getCompFrac(triplet.comps[j], triplet.cells[j], "pec");
        if (pecFrac == 0)
          use = false;
      }
    }

    if (use)
      triplets.push_back(triplet);
  }
}
#endif

// careful when averaging the normal for multiple mu objects
template<size_t DIM, typename Scalar>
MxDimVector<double, DIM> MxYeeFitMu<DIM, Scalar>::getNormal(
typename MxYeeFitMu<DIM, Scalar>::Stencil const & stencil) const {
  std::vector<MxDimVector<double, DIM> > norms;
  int count = 0;

  MxMu<DIM> const * muObj;
  double lfrac, afrac;
  bool cut;
  for (size_t i = 0; i < sim->numMus(); ++i) {
    muObj = sim->getMu(i);
    cut = false;

    for (size_t j = 0; j < stencil.comps.size(); ++j) {
      lfrac = hfield->getCompFrac(stencil.comps[j], stencil.cells[j],
        muObj->getName());
      afrac = bfield->getCompFrac(stencil.comps[j], stencil.cells[j],
        muObj->getName());
      if ((lfrac != 0 and lfrac != 1) or (afrac != 0 and afrac != 1))
        cut = true;
    }

    if (cut)
      norms.push_back(muObj->getShape()->normal(hfield->getCompCoord(
        stencil.comps[0], stencil.cells[0])));
  }

  // Relevant normals gathered. Do the averaging.
  MxDimVector<double, DIM> norm(0);
  if (norms.size() > 0) {
    // Adjust the sign of each norm (except for the first) such that each
    // one's dot product with the first is nonnegative
    std::vector<double> signs(norms.size(), 1);
    double dotProd;
    for (size_t i = 1; i < signs.size(); ++i) {
      dotProd = norms[0].dot(norms[i]);
      if (dotProd < 0)
        signs[i] = -1;
    }

    // add
    for (size_t i = 0; i < signs.size(); ++i)
      norm += signs[i] * norms[i];

    return norm / norm.norm();
  }
  else {
    // if there are no cuts in the current stencil, any normal will do.
    norm[0] = 1;
    return norm;
  }

}

// multi tuple update. This function is meant to handle several different
// mu objects with possibly complex mu constants.
template<size_t DIM, typename Scalar>
MxDimMatrix<MxComplex, DIM> MxYeeFitMu<DIM, Scalar>::tupleUpdate(
typename MxYeeFitMu<DIM, Scalar>::Tuple const & tuple,
MxDimVector<double, DIM> normal) {

  // need the normal vector as a complex vector
  MxDimVector<MxComplex, DIM> normc;
  for (size_t i = 0; i < DIM; ++i)
    normc[i] = MxComplex(normal[i], 0.0);

  // complex cartesian projection matrices
  std::vector<MxDimMatrix<MxComplex, DIM> > cartProjs(DIM,
    MxDimMatrix<MxComplex, DIM>(0));
  for (size_t i = 0; i < DIM; ++i)
    cartProjs[i](i, i) = ScalarTraits<MxComplex>::one();

  MxDimMatrix<MxComplex, DIM> muDim, gamma, pi,
      aveGamma(ScalarTraits<MxComplex>::zero()),
      avePi(ScalarTraits<MxComplex>::zero());
  MxDimMatrix<MxComplex, DIM> nn(normc, normc),
      eye(MxDimMatrix<MxComplex, DIM>::I());

  double lfrac, afrac;
  MxDimVector<double, DIM> lfracSums(0), afracSums(0);
  MxDimVector<MxComplex, DIM> lfracs, afracs;
  MxMu<DIM> const * muObj;
  MxDimMatrix<MxComplex, 3> mu;
  size_t comp;
  for (size_t i = 0; i < sim->numMus(); ++i) {
    muObj = sim->getMu(i);

    for (size_t j = 0; j < DIM; ++j) {
      comp = tuple.comps[j];

      lfrac = hfield->getCompFrac(comp, tuple.cells[j], muObj->getName());
      lfracSums[comp] += lfrac;
      lfracs[comp] = MxComplex(lfrac, 0);

      afrac = bfield->getCompFrac(comp, tuple.cells[j], muObj->getName());
      afracSums[comp] += afrac;
      afracs[comp] = MxComplex(afrac, 0);
    }

    mu = muObj->getMu();
    for (size_t j = 0; j < DIM; ++j)
      for (size_t k = 0; k < DIM; ++k)
        muDim(j, k) = mu(j, k);

    gamma = eye + nn * (eye - muDim) / normc.dot(muDim * normc);
    pi = muDim * gamma;

    for (size_t j = 0; j < DIM; ++j) {
      aveGamma += cartProjs[j] * lfracs[j] * gamma;
      avePi += cartProjs[j] * afracs[j] * pi;
    }
  }

  // now add background mu
  muObj = sim->getBackgroundMu();

  for (size_t j = 0; j < DIM; ++j) {
    lfrac = 1.0 - lfracSums[j];
    lfracs[j] = MxComplex(lfrac, 0);

    afrac = 1.0 - afracSums[j];
    afracs[j] = MxComplex(afrac, 0);
  }

  mu = muObj->getMu();
  for (size_t j = 0; j < DIM; ++j)
    for (size_t k = 0; k < DIM; ++k)
      muDim(j, k) = mu(j, k);

  gamma = eye + nn * (eye - muDim) / normc.dot(muDim * normc);
  pi = muDim * gamma;

  for (size_t j = 0; j < DIM; ++j) {
    aveGamma += cartProjs[j] * lfracs[j] * gamma;
    avePi += cartProjs[j] * afracs[j] * pi;
  }

#if 0 
  std::cout << "lfracs: "; lfracs.print();
  std::cout << "afracs: "; afracs.print();
  std::cout << "normal: "; normal.print();
  std::cout << "update matrix\n";
  (aveGamma * avePi.inv()).print();
#endif

  // SWITCH
  if (mInvert)
    return aveGamma * avePi.inv(); // this is for mu^{-1}
  else
    return avePi * aveGamma.inv(); // this is for mu
}



template<size_t DIM, typename Scalar>
void MxYeeFitMu<DIM, Scalar>::setMatrix() {

  MxGridFieldIter<DIM> fieldIter(hfield);

  MxIndex row;
  size_t comp0;
  MxDimVector<int, DIM> cell;
  MxDimVector<double, DIM> coord;

  double lfrac0, afrac0;

  // Mu stuff
  bool hasMu = sim->hasMu();
  MxMu<DIM> const * muObj(0);
  bool inMu;
  bool muIsDiag;

  // PML stuff
  bool hasPml = sim->hasPML();
  MxPML<DIM> const * pml(0);
  double fPml;
  bool inPml;
  bool pmlIsDiag;
  MxDimMatrix<MxComplex, 3> pmlTransf(MxDimMatrix<MxComplex, 3>::I());
  MxDimMatrix<MxComplex, DIM> pmlTransfDim(MxDimMatrix<MxComplex, DIM>::I());
  MxComplex detPmlTransf(1.0, 0.0);

  // for the complicated updates...
  size_t stencilSize = 1 + (DIM - 1) * 4;  // =5 in 2D, =9 in 3D
  typename MxYeeFitMu<DIM, Scalar>::Stencil stencil;
  std::vector<typename MxYeeFitMu<DIM, Scalar>::Tuple> tuples;
  size_t numTuples;
  double numTuplesDbl;
  MxDimVector<double, DIM> n; // normal direction for mu update
  MxDimMatrix<MxComplex, DIM> update; // update matrix for one tuple
  size_t stencilInd;
  size_t comp;

  // for diagonal updates...
  MxComplex factor;
  MxComplex cval;
  Scalar val;

  // debug
  int countDiag = 0;
  int countFull = 0;

  // TE: Bz, Ex, Ey
  for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
    comp0 = fieldIter.getComp(); // Ex : 0; Ey : 1
    cell = fieldIter.getCell();
    row = fieldIter.getGlobCompIndx();
    coord = fieldIter.getCoord();


    // check for trivial mu updates. Assume nontrivial at first
    inMu = false;
    muIsDiag = false;
    if (hasMu) {
      // loop through all mu objects in simulation
      for (size_t i = 0; i < sim->numMus(); ++i) {
        muObj = sim->getMu(i);
        lfrac0 = hfield->getCompFrac(comp0, cell, muObj->getName());
        afrac0 = bfield->getCompFrac(comp0, cell, muObj->getName()); 

        // first case: we are completely inside the current
        // mu and the muilon tensor is diagonal. This causes
        // the loop over mus to end, since we assume that
        // there are no overlapping mu objects.
        if (lfrac0 == 1 and afrac0 == 1) {
          inMu = true;
          muIsDiag = muObj->isMuDiag();
          break;
        }
        // second case: outside mu, do nothing
        else if (lfrac0 == 0 and afrac0 == 0)
          continue;
        // third case: lfrac or afrac is cut. Call this in the
        // mu and assume muilon is NOT diagonal
        else {
          inMu = true;
          break;
        }
      }
    }
    // get background mu if component wasn't inside any
    // object
    if (not inMu) {
      muObj = sim->getBackgroundMu();
      muIsDiag = muObj->isMuDiag();
    }


    // check for intersection with pmls
    // by default, assume pml is diagonal
    pmlIsDiag = true;
    if (sim->hasPML()) {
      // reset pmlTransf
      pmlTransf = MxDimMatrix<MxComplex, 3>::I();

      // loop through all pml objects in simulation
      for (size_t i = 0; i < sim->numPMLs(); ++i) {
        pml = sim->getPML(i);

        fPml = pml->getShape()->func(coord);
        if (fPml > 0) {
          if (not pml->isDiag())
            pmlIsDiag = false;

          if (mInvert)
            pmlTransf = pmlTransf * pml->getInvS(coord);
          else
            pmlTransf = pmlTransf * pml->getS(coord);
        }
      }
      detPmlTransf = det(pmlTransf);

      for (size_t j = 0; j < DIM; ++j)
        for (size_t k = 0; k < DIM; ++k)
          pmlTransfDim(j, k) = pmlTransf(j, k);
      //std::cout << "coord: " << coord;
      //std::cout << pmlTransf;
      //std::cout << detPmlTransf << "\n";
    }


    if (muIsDiag and pmlIsDiag) {
      //countDiag++;
      factor = bfield->getCompFactor(comp0, cell);
      if (mInvert)
        cval = factor / muObj->getMu()(comp0, comp0);
      else
        cval = factor * muObj->getMu()(comp0, comp0);
      cval *= pmlTransf(comp0, comp0) * pmlTransf(comp0, comp0) / detPmlTransf;

      MxUtil::convertScalar(cval, val);
      MxCrsMatrix<Scalar>::insertRowValues(row, 1, &row, &val);
    }
    else {
      //countFull++;
      stencil = getStencil(comp0, cell);

      // averaged normal direction for mus
      n = getNormal(stencil);
      //std::cout << n;

      // pick out tuples from stencil
      getTuples(stencil, tuples);
      numTuples = tuples.size();
      numTuplesDbl = double(numTuples);

      // loop through the tuples, calculating the mu update for each
      // and adding the result to the update stencil values.
      for (size_t i = 0; i < numTuples; ++i) {
        update = tupleUpdate(tuples[i], n);
        // PML transformation of update...
        update = (pmlTransfDim * update * pmlTransfDim) / detPmlTransf; 
      
        for (size_t j = 0; j < DIM; ++j) {
          stencilInd = tuples[i].stencilInds[j];
          comp = tuples[i].comps[j];
          MxUtil::convertScalar(update(comp0, comp) / numTuplesDbl, val);
          stencil.values[stencilInd] += val;
        }
      }

      // enforce boundary conditions (i.e. zero for PEC bcs, phase shifts for
      // periodic simulation...)
      for (size_t i = 0; i < stencilSize; ++i)
        stencil.values[i] *= MxUtil::convertScalar(stencil.bcfactors[i], val);

      //MxUtil::printStdVector(stencil.rvalues);
      //MxUtil::printStdVector(stencil.columns);

      // put 'em on in there
      MxCrsMatrix<Scalar>::insertRowValues(
        row, stencil.columns, stencil.values);
    }

  }

  //std::cout << "Diagonal updates: " << countDiag << "\n";
  //std::cout << "Full updates: " << countFull << "\n";
  MxCrsMatrix<Scalar>::fillComplete(bfield->getMap(), hfield->getMap());

}



template<size_t DIM, typename Scalar>
void MxYeeFitMu<DIM, Scalar>::setMatrix2dTE() {

  MxDimVector<int, DIM> cell;
  MxDimVector<double, DIM> coord;

  MxIndex row;
  Scalar val;
  double afrac;

  double sumAfracs;
  MxComplex aveMu, cval;

  MxDimMatrix<MxComplex, 3> s;

  MxMu<DIM> const * muObj;
  MxPML<DIM> const * pml;

  // TE: Hz, Ex, Ey
  MxGridFieldIter<DIM> fieldIter(hfield);
  for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
    cell = fieldIter.getCell();
    row = hfield->globCompIndx(0, cell);
    coord = fieldIter.getCoord();

    aveMu = 0;
    sumAfracs = 0;
    for (size_t i = 0; i < sim->numMus(); ++i) {
      muObj = sim->getMu(i);
      afrac = bfield->getCompFrac(0, cell, muObj->getName()); 

      aveMu += afrac * muObj->getMu()(2, 2);
      sumAfracs += afrac;
    }
    aveMu += (1.0 - sumAfracs) * sim->getBackgroundMu()->getMu()(2, 2);

    // now do pml transformation
    if (sim->hasPML()) {
      // reset pml transformation
      s = MxDimMatrix<MxComplex, 3>::I();

      // loop through all pml objects in simulation
      for (size_t i = 0; i < sim->numPMLs(); ++i) {
        pml = sim->getPML(i);

        if (pml->getShape()->func(coord) > 0)
          s = s * pml->getS(coord);
      }
      //std::cout << "coord: " << coord;
      //std::cout << pmlTransf;
      //std::cout << detPmlTransf << "\n";
      aveMu *= s(2,2)*s(2,2) / det(s);
    }

    if (mInvert)
      cval = 1.0 / aveMu;
    else
      cval = aveMu;

    MxUtil::convertScalar(cval, val);
    MxCrsMatrix<Scalar>::insertRowValues(row, 1, &row, &val);
  }

  MxCrsMatrix<Scalar>::fillComplete(bfield->getMap(), hfield->getMap());

}

// returns pointer to matrix. Receiver is responsible for deleting it.
template<size_t DIM, typename Scalar>
RCP<MxCrsMatrix<Scalar> > MxYeeFitMu<DIM, Scalar>::cellAveInvMuOperator() const {
  size_t numEComps = hfield->getNumComps();

  std::vector<MxComplex> scalarMus;
  MxDimMatrix<MxComplex, 3> mu;
  for (size_t i = 0; i < sim->numMus(); ++i) {
    mu = sim->getMu(i)->getMu();
    if (DIM == 2 and numEComps == 1)
      scalarMus.push_back(1.0 / mu(2, 2));
    else if (DIM == 2 and numEComps == 2)
      scalarMus.push_back(2.0 / (mu(0, 0) + mu(1, 1)));
    else
      scalarMus.push_back(3.0 / mu.trace());
  }
  MxComplex bgScalarMu;
  mu = sim->getBackgroundMu()->getMu();
  if (DIM == 2 and numEComps == 1)
    bgScalarMu = 1.0 / mu(2, 2);
  else if (DIM == 2 and numEComps == 2)
    bgScalarMu = 2.0 / (mu(0, 0) + mu(1, 1));
  else
    bgScalarMu = 3.0 / mu.trace();


  MxPML<DIM> const * pml(0);
  MxDimMatrix<MxComplex, DIM> pmlTransf(MxDimMatrix<MxComplex, DIM>::I());

  const MxGridField<DIM> & psifield = sim->getField("psifield");
  RCP<MxCrsMatrix<Scalar> > res;
  res = rcp(new MxCrsMatrix<Scalar>(psifield.getMap(), psifield.getMap()));

  MxGridFieldIter<DIM> iter(&psifield);
  double frac, sumFracs;
  MxIndex row;
  size_t comp;
  MxDimVector<int, DIM> cell;
  MxDimVector<double, DIM> coord;
  MxComplex aveInvMu;
  Scalar val;

  for (iter.begin(); !iter.atEnd(); iter.bump()) {
    row = iter.getGlobCompIndx();
    comp = iter.getComp();
    cell = iter.getCell();

    sumFracs = 0;
    aveInvMu = 0;
    for (size_t i = 0; i < sim->numMus(); ++i) {
      frac = psifield.getCompFrac(comp, cell,
        sim->getMu(i)->getName());
      sumFracs += frac;

      aveInvMu += frac * scalarMus[i];
    }
    aveInvMu += (1.0 - sumFracs) * bgScalarMu;

#if 0
    // now get pml 'volume'
    if (sim->hasPML()) {
      // reset invS
      pmlTransf = MxDimMatrix<MxComplex, DIM>::I();

      // loop through all pml objects in simulation
      for (size_t i = 0; i < sim->numPMLs(); ++i) {
        pml = sim->getPML(i);
        if (pml->getShape()->func(coord) > 0)
          pmlTransf = pmlTransf * pml->getInvS(coord);
      }
      aveInvMu *= det(pmlTransf);
    }
#endif

    MxUtil::convertScalar(aveInvMu, val);
    res->insertRowValues(row, 1, &row, &val);
  }

  res->fillComplete(psifield.getMap(), psifield.getMap());
  return res;
}


template class MxYeeFitMu<1, double>;
template class MxYeeFitMu<2, double>;
template class MxYeeFitMu<3, double>;
template class MxYeeFitMu<1, MxComplex>;
template class MxYeeFitMu<2, MxComplex>;
template class MxYeeFitMu<3, MxComplex>;
