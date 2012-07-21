#include "MxYeeFitInvEps.h"

#include <limits>

#include "MxDimVector.hpp"
#include "MxGridFieldIter.hpp"
#include "MxIO.h"



template<size_t DIM, typename Scalar>
MxYeeFitInvEps<DIM, Scalar>::MxYeeFitInvEps(RCP<MxEMSim<DIM> > theSim,
bool invert) :
//MxCrsMatrix<Scalar>(theSim->getField("efield").getMap(),
//  theSim->getField("efield").getMap()),
MxCrsMatrix<Scalar>(theSim->getField("efield").getMap()),
sim(theSim), efield(&theSim->getField("efield")),
dfield(&theSim->getField("dfield")), mInvert(invert),
grid(&theSim->getGrid()) {

  if (efield->getNumComps() == 1)
    setMatrix2dTM();
  else
    setMatrix();

  //MxPML<DIM>::alphaMap(theSim, efield);
  saveEps();
}




template<size_t DIM, typename Scalar>
typename MxYeeFitInvEps<DIM, Scalar>::Stencil
MxYeeFitInvEps<DIM, Scalar>::getStencil(size_t comp0,
MxDimVector<int, DIM> cell) const {
  size_t comp1 = (comp0 + 1) % efield->getNumComps();
  size_t comp2 = (comp1 + 1) % efield->getNumComps();

  size_t stenSize = 1 + 4*(DIM - 1); // 5 in 2D; 9 in 3D

  typename MxYeeFitInvEps<DIM, Scalar>::Stencil stencil;
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
  stencil.cells[1][comp1]--;
  stencil.cells[3][comp0]++; stencil.cells[3][comp1]--;
  stencil.cells[4][comp0]++;
  if (DIM == 3) {
    stencil.cells[5][comp2]--;
    stencil.cells[7][comp0]++; stencil.cells[7][comp2]--;
    stencil.cells[8][comp0]++;
  }

  stencil.columns = std::vector<MxIndex>(stenSize);
  for (size_t i = 0; i < stenSize; ++i)
    stencil.columns[i] =
      dfield->globCompIndx(stencil.comps[i], stencil.cells[i]);

  stencil.bcfactors = std::vector<MxComplex>(stenSize);
  for (size_t i = 0; i < stenSize; ++i)
    stencil.bcfactors[i] =
      dfield->getCompFactor(stencil.comps[i], stencil.cells[i]);

  return stencil;
}

// obsolete due to above generic version
#if 0
typedef MxYeeFitInvEps<2>::Stencil Stencil2D;
typedef MxYeeFitInvEps<3>::Stencil Stencil3D;

template<>
MxYeeFitInvEps<1>::Stencil MxYeeFitInvEps<1>::getStencil(int comp0,
MxDimVector<int, 1> cell) const {
  return MxYeeFitInvEps<1>::Stencil();
}

template<>
Stencil3D MxYeeFitInvEps<3>::getStencil(int comp0,
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
  stencil.cells[1][comp1]--;
  stencil.cells[3][comp0]++; stencil.cells[3][comp1]--;
  stencil.cells[4][comp0]++;
  stencil.cells[5][comp2]--;
  stencil.cells[7][comp0]++; stencil.cells[7][comp2]--;
  stencil.cells[8][comp0]++;

  stencil.columns = std::vector<int>(9);
  for (size_t i = 0; i < 9; ++i)
    stencil.columns[i] = dfield->globCompIndx(stencil.comps[i], stencil.cells[i]);

  stencil.bcfactors = std::vector<MxComplex>(9);
  for (size_t i = 0; i < 9; ++i)
    stencil.bcfactors[i] = dfield->getCompFactor(stencil.comps[i], stencil.cells[i]);

  return stencil;
}
#endif

template<size_t DIM, typename Scalar>
void MxYeeFitInvEps<DIM, Scalar>::getTuples(
typename MxYeeFitInvEps<DIM, Scalar>::Stencil const & stencil,
std::vector<typename MxYeeFitInvEps<DIM, Scalar>::Tuple> & tuples) const {
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
  
  typename MxYeeFitInvEps<DIM, Scalar>::Tuple tuple;

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
        pecFrac = efield->getCompFrac(tuple.comps[j], tuple.cells[j], "pec");
        if (pecFrac == 0)
          use = false;
      }
    }

    if (use)
      tuples.push_back(tuple);
  }
}

// obsolete! yey!
#if 0
typedef MxYeeFitInvEps<2>::Tuple Doublet;
typedef MxYeeFitInvEps<3>::Tuple Triplet;

template<>
void MxYeeFitInvEps<1>::getTuples(MxYeeFitInvEps<1>::Stencil const & stencil,
std::vector<MxYeeFitInvEps<1>::Tuple> & singlets) const {}

template<>
void MxYeeFitInvEps<2>::getTuples(Stencil2D const & stencil,
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
        pecFrac = efield->getCompFrac(doublet.comps[j], doublet.cells[j], "pec");
        if (pecFrac == 0)
          use = false;
      }
    }

    if (use)
      doublets.push_back(doublet);
  }
}

template<>
void MxYeeFitInvEps<3>::getTuples(Stencil3D const & stencil,
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
        pecFrac = efield->getCompFrac(triplet.comps[j], triplet.cells[j], "pec");
        if (pecFrac == 0)
          use = false;
      }
    }

    if (use)
      triplets.push_back(triplet);
  }
}
#endif


// careful when averaging the normal for multiple dielectric objects
template<size_t DIM, typename Scalar>
MxDimVector<double, DIM> MxYeeFitInvEps<DIM, Scalar>::getNormal(
MxYeeFitInvEps<DIM, Scalar>::Stencil const & stencil) const {
  std::vector<MxDimVector<double, DIM> > norms;
  int count = 0;

  MxDielectric<DIM> const * diel;
  double lfrac, afrac;
  bool cut;
  for (size_t i = 0; i < sim->numDielectrics(); ++i) {
    diel = sim->getDielectric(i);
    cut = false;

    for (size_t j = 0; j < stencil.comps.size(); ++j) {
      lfrac = efield->getCompFrac(stencil.comps[j], stencil.cells[j],
        diel->getName());
      afrac = dfield->getCompFrac(stencil.comps[j], stencil.cells[j],
        diel->getName());
      if ((lfrac != 0 and lfrac != 1) or (afrac != 0 and afrac != 1))
        cut = true;
    }

    if (cut)
      norms.push_back(diel->getShape()->normal(efield->getCompCoord(
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
// dielectric objects with possibly complex dielectric constants.
template<size_t DIM, typename Scalar>
MxDimMatrix<MxComplex, DIM> MxYeeFitInvEps<DIM, Scalar>::tupleUpdate(
MxYeeFitInvEps<DIM, Scalar>::Tuple const & tuple,
MxDimVector<double, DIM> normal) {

  // need the normal vector as a complex vector
  MxDimVector<MxComplex, DIM> normc;
  for (size_t i = 0; i < DIM; ++i)
    normc[i] = MxComplex(normal[i], 0.0);

  // complex cartesian projection matrices
  std::vector<MxDimMatrix<MxComplex, DIM> > cartProjs(
    DIM, MxDimMatrix<MxComplex, DIM>(0));
  for (size_t i = 0; i < DIM; ++i)
    cartProjs[i](i, i) = ScalarTraits<MxComplex>::one();

  MxDimMatrix<MxComplex, DIM> epsDim, gamma, pi,
      aveGamma(ScalarTraits<MxComplex>::zero()),
      avePi(ScalarTraits<MxComplex>::zero());
  MxDimMatrix<MxComplex, DIM> nn(normc, normc),
      eye(MxDimMatrix<MxComplex, DIM>::I());

  double lfrac, afrac;
  MxDimVector<double, DIM> lfracSums(0), afracSums(0);
  MxDimVector<MxComplex, DIM> lfracs, afracs;
  MxDielectric<DIM> const * diel;
  MxDimMatrix<MxComplex, 3> eps;
  size_t comp;
  for (size_t i = 0; i < sim->numDielectrics(); ++i) {
    diel = sim->getDielectric(i);

    for (size_t j = 0; j < DIM; ++j) {
      comp = tuple.comps[j];

      lfrac = efield->getCompFrac(comp, tuple.cells[j], diel->getName());
      lfracSums[comp] += lfrac;
      lfracs[comp] = MxComplex(lfrac, 0);

      afrac = dfield->getCompFrac(comp, tuple.cells[j], diel->getName());
      afracSums[comp] += afrac;
      afracs[comp] = MxComplex(afrac, 0);
    }

    eps = diel->getEps();
    for (size_t j = 0; j < DIM; ++j)
      for (size_t k = 0; k < DIM; ++k)
        epsDim(j, k) = eps(j, k);

    gamma = eye + nn * (eye - epsDim) / normc.dot(epsDim * normc);
    pi = epsDim * gamma;

    for (size_t j = 0; j < DIM; ++j) {
      aveGamma += cartProjs[j] * lfracs[j] * gamma;
      avePi += cartProjs[j] * afracs[j] * pi;
    }
  }

  // now add background eps
  diel = sim->getBackgroundDielectric();

  for (size_t j = 0; j < DIM; ++j) {
    lfrac = 1.0 - lfracSums[j];
    lfracs[j] = MxComplex(lfrac, 0);

    afrac = 1.0 - afracSums[j];
    afracs[j] = MxComplex(afrac, 0);
  }

  eps = diel->getEps();
  for (size_t j = 0; j < DIM; ++j)
    for (size_t k = 0; k < DIM; ++k)
      epsDim(j, k) = eps(j, k);

  gamma = eye + nn * (eye - epsDim) / normc.dot(epsDim * normc);
  pi = epsDim * gamma;

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
  return aveGamma * avePi.inv();
}



template<size_t DIM, typename Scalar>
void MxYeeFitInvEps<DIM, Scalar>::setMatrix() {

  MxGridFieldIter<DIM> fieldIter(efield);

  MxIndex row;
  size_t comp0;
  MxDimVector<int, DIM> cell;
  MxDimVector<double, DIM> coord;

  double lfrac0, afrac0;

  // Dielectric stuff
  bool hasDiel = sim->hasDielectric();
  MxDielectric<DIM> const * diel(0);
  bool inDiel;
  bool epsIsDiag;

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
  typename MxYeeFitInvEps<DIM, Scalar>::Stencil stencil;
  std::vector<typename MxYeeFitInvEps<DIM, Scalar>::Tuple> tuples;
  size_t numTuples;
  double numTuplesDbl;
  MxDimVector<double, DIM> n; // normal direction for dielectric update
  MxDimMatrix<MxComplex, DIM> invEps; // update matrix for one tuple
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


    // check for trivial dielectric updates. Assume nontrivial at first
    inDiel = false;
    epsIsDiag = false;
    if (hasDiel) {
      // loop through all dielectric objects in simulation
      for (size_t i = 0; i < sim->numDielectrics(); ++i) {
        diel = sim->getDielectric(i);
        lfrac0 = efield->getCompFrac(comp0, cell, diel->getName());
        afrac0 = dfield->getCompFrac(comp0, cell, diel->getName()); 

        // first case: we are completely inside the current
        // dielectric and the epsilon tensor is diagonal. This causes
        // the loop over dielectrics to end, since we assume that
        // there are no overlapping dielectric objects.
        if (lfrac0 == 1 and afrac0 == 1) {
          inDiel = true;
          epsIsDiag = diel->isEpsDiag();
          break;
        }
        // second case: outside dielectric, do nothing
        else if (lfrac0 == 0 and afrac0 == 0)
          continue;
        // third case: lfrac or afrac is cut. Call this in the
        // dielectric and assume epsilon is NOT diagonal
        else {
          inDiel = true;
          break;
        }
      }
    }
    // get background dielectric if component wasn't inside any
    // object
    if (not inDiel) {
      diel = sim->getBackgroundDielectric();
      epsIsDiag = diel->isEpsDiag();
    }


    // check for intersection with pmls
    // by default, assume pml is diagonal
    pmlIsDiag = true;
    if (sim->hasPML()) {
      // reset invS
      pmlTransf = MxDimMatrix<MxComplex, 3>::I();

      // loop through all pml objects in simulation
      for (size_t i = 0; i < sim->numPMLs(); ++i) {
        pml = sim->getPML(i);

        fPml = pml->getShape()->func(coord);
        if (fPml > 0) {
          if (not pml->isDiag())
            pmlIsDiag = false;
          pmlTransf = pmlTransf * pml->getInvS(coord);
        }
      }
      detPmlTransf = det(pmlTransf);

      for (size_t j = 0; j < DIM; ++j)
        for (size_t k = 0; k < DIM; ++k)
          pmlTransfDim(j, k) = pmlTransf(j, k);

      //std::cout << pmlTransfDim;
    }


    if (epsIsDiag and pmlIsDiag) {
      //countDiag++;
      factor = dfield->getCompFactor(comp0, cell);
      cval = factor * pmlTransf(comp0, comp0) * pmlTransf(comp0, comp0) / 
        detPmlTransf / diel->getEps()(comp0, comp0);
      MxUtil::convertScalar(cval, val);
      MxCrsMatrix<Scalar>::insertRowValues(row, 1, &row, &val);
    }
    else {
      //countFull++;
      stencil = getStencil(comp0, cell);

      // averaged normal direction for dielectrics
      n = getNormal(stencil);
      //std::cout << n;

      // pick out tuples from stencil
      getTuples(stencil, tuples);
      numTuples = tuples.size();
      numTuplesDbl = double(numTuples);

      // loop through the tuples, calculating the dielectric update for each
      // and adding the result to the update stencil values.
      for (size_t i = 0; i < numTuples; ++i) {
        invEps = tupleUpdate(tuples[i], n); // dielectric part of invEps
        invEps = (pmlTransfDim * invEps * pmlTransfDim) / detPmlTransf; // PML transformation
      
        for (size_t j = 0; j < DIM; ++j) {
          stencilInd = tuples[i].stencilInds[j];
          comp = tuples[i].comps[j];
          MxUtil::convertScalar(invEps(comp0, comp) / numTuplesDbl, val);
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
  MxCrsMatrix<Scalar>::fillComplete(efield->getMap(), dfield->getMap());

}

//template<size_t DIM>
//void MxYeeFitInvEps<DIM>::build() {
//  for (int i = 0; i < num
//
//}

template<size_t DIM, typename Scalar>
void MxYeeFitInvEps<DIM, Scalar>::saveEps() const {
  size_t numComps = efield->getNumComps();

  std::cout << "create eps vec\n";
  MxMultiVector<MxComplex> res(efield->getMap(), 1);
  std::cout << "done\n";

  MxGrid<DIM> const & grid = sim->getGrid();

  bool hit;
  size_t comp;
  MxIndex row;
  MxDimVector<double, DIM> coord;
  MxDimMatrix<MxComplex, 3> eps;
  MxDielectric<DIM> const * diel;
  
  MxGridFieldIter<DIM> fieldIter(efield);
  for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
    if (DIM == 2 and numComps == 1)
      comp = 2;
    else
      comp = fieldIter.getComp();

    coord = fieldIter.getCoord();
    row = fieldIter.getGlobCompIndx();
    std::cout << "row=" << row << "\n";

    hit = false;
    for (size_t i = 0; i < sim->numDielectrics(); ++i) {
      diel = sim->getDielectric(i);
      if (diel->getShape()->func(coord) > 0) {
        hit = true;
        res.replaceGlobalValue(row, 0, diel->getEps()(comp, comp));
      }
    }
    if (hit == false) {
      eps = sim->getBackgroundDielectric()->getEps();
      res.replaceGlobalValue(row, 0, eps(comp, comp));
    }
  }

  MxIO<DIM, MxComplex> io(sim->getGrid().getComm());
  std::cout << "Saving eps\n";
  io.save(res, *efield, "epsilon");

}

template<size_t DIM, typename Scalar>
void MxYeeFitInvEps<DIM, Scalar>::setMatrix2dTM() {

  MxDimVector<int, DIM> cell;
  MxDimVector<double, DIM> coord;

  MxIndex row;
  Scalar val;
  double afrac;

  double sumAfracs;
  MxComplex aveEps, cval;

  MxDimMatrix<MxComplex, 3> s;

  MxDielectric<DIM> const * epsObj;
  MxPML<DIM> const * pml;

  // TE: Ez, Hx, Hy
  MxGridFieldIter<DIM> fieldIter(efield);
  for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
    cell = fieldIter.getCell();
    row = efield->globCompIndx(0, cell);
    coord = fieldIter.getCoord();

    aveEps = 0;
    sumAfracs = 0;
    for (size_t i = 0; i < sim->numDielectrics(); ++i) {
      epsObj = sim->getDielectric(i);
      afrac = dfield->getCompFrac(0, cell, epsObj->getName()); 

      aveEps += afrac * epsObj->getEps()(2, 2);
      sumAfracs += afrac;
    }
    aveEps += (1.0 - sumAfracs) * sim->getBackgroundDielectric()->getEps()(2, 2);

    // now do pml transformation
    if (sim->hasPML()) {
      // reset pml transformation
      s = MxDimMatrix<MxComplex, 3>::I();

      // loop through all pml objects in siepslation
      for (size_t i = 0; i < sim->numPMLs(); ++i) {
        pml = sim->getPML(i);

        if (pml->getShape()->func(coord) > 0)
          s = s * pml->getS(coord);
      }
      //std::cout << "coord: " << coord;
      //std::cout << pmlTransf;
      //std::cout << detPmlTransf << "\n";
      aveEps *= s(2,2)*s(2,2) / det(s);
    }

    if (mInvert)
      cval = 1.0 / aveEps;
    else
      cval = aveEps;

    MxUtil::convertScalar(cval, val);
    MxCrsMatrix<Scalar>::insertRowValues(row, 1, &row, &val);
  }

  MxCrsMatrix<Scalar>::fillComplete(dfield->getMap(), efield->getMap());
}


// returns pointer to matrix. Receiver is responsible for deleting it.
template<size_t DIM, typename Scalar>
RCP<MxCrsMatrix<Scalar> > MxYeeFitInvEps<DIM, Scalar>::cellAveInvEpsOperator() const {
  size_t numEComps = efield->getNumComps();

  // get scalar epsilons
  std::vector<MxComplex> scalarEpss;
  MxDimMatrix<MxComplex, 3> eps;
  for (size_t i = 0; i < sim->numDielectrics(); ++i) {
    eps = sim->getDielectric(i)->getEps();
    if (DIM == 2 and numEComps == 1)
      scalarEpss.push_back(1.0 / eps(2, 2));
    else if (DIM == 2 and numEComps == 2)
      scalarEpss.push_back(2.0 / (eps(0, 0) + eps(1, 1)));
    else
      scalarEpss.push_back(3.0 / eps.trace());
  }
  MxComplex bgScalarEps;
  eps = sim->getBackgroundDielectric()->getEps();
  if (DIM == 2 and numEComps == 1)
    bgScalarEps = 1.0 / eps(2, 2);
  else if (DIM == 2 and numEComps == 2)
    bgScalarEps = 2.0 / (eps(0, 0) + eps(1, 1));
  else
    bgScalarEps = 3.0 / eps.trace();


  MxPML<DIM> const * pml(0);
  MxDimMatrix<MxComplex, 3> pmlTransf(MxDimMatrix<MxComplex, 3>::I());

  const MxGridField<DIM> & psifield = sim->getField("psifield");
  RCP<MxCrsMatrix<Scalar> > res;
  res = rcp(new MxCrsMatrix<Scalar>(psifield.getMap(), psifield.getMap()));

  MxGridFieldIter<DIM> iter(&psifield);
  double frac, sumFracs;
  MxIndex row;
  size_t comp;
  MxDimVector<int, DIM> cell;
  MxDimVector<double, DIM> coord;
  MxComplex aveInvEps;
  Scalar val;

  for (iter.begin(); !iter.atEnd(); iter.bump()) {
    row = iter.getGlobCompIndx();
    comp = iter.getComp();
    cell = iter.getCell();
    coord = iter.getCoord();

    sumFracs = 0;
    aveInvEps = 0;
    for (size_t i = 0; i < sim->numDielectrics(); ++i) {
      frac = psifield.getCompFrac(comp, cell,
        sim->getDielectric(i)->getName());
      sumFracs += frac;

      aveInvEps += frac * scalarEpss[i];
    }
    aveInvEps += (1.0 - sumFracs) * bgScalarEps;

    // now get pml 'volume'
    if (sim->hasPML()) {
      // reset invS
      pmlTransf = MxDimMatrix<MxComplex, 3>::I();

      // loop through all pml objects in simulation
      for (size_t i = 0; i < sim->numPMLs(); ++i) {
        pml = sim->getPML(i);
        if (pml->getShape()->func(coord) > 0)
          pmlTransf = pmlTransf * pml->getS(coord);
      }
      aveInvEps *= det(pmlTransf);
    }

    MxUtil::convertScalar(aveInvEps, val);
    res->insertRowValues(row, 1, &row, &val);
  }

  res->fillComplete(psifield.getMap(), psifield.getMap());
  return res;
}


template class MxYeeFitInvEps<1, double>;
template class MxYeeFitInvEps<2, double>;
template class MxYeeFitInvEps<3, double>;
template class MxYeeFitInvEps<1, MxComplex>;
template class MxYeeFitInvEps<2, MxComplex>;
template class MxYeeFitInvEps<3, MxComplex>;
