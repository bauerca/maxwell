
#include "MxYeeMitInvEps.h"

#include <limits>

#include "MxDimVector.hpp"
#include "MxGridFieldIter.hpp"
#include "MxCartRect.hpp"
#include "MxCartBox.hpp"



template<size_t DIM>
MxYeeMitInvEps<DIM>::MxYeeMitInvEps(const MxEMSim<DIM> * theSim, bool imag) :
Epetra_CrsMatrix(Copy, theSim->getField("efield").getMap(), 0),
sim(theSim), efield(&theSim->getField("efield")), dfield(&theSim->getField("dfield")),
grid(&theSim->getGrid()) {
  this->imaginary = imag;
  if (efield->getNumComps() == 1)
    setMatrix2dTM();
  else {
    std::cout << "filling inv eps's\n";
    fillAveEps();
    std::cout << "building matrix\n";
    setMatrix2dTE3d();
  }
}


template<size_t DIM>
bool MxYeeMitInvEps<DIM>::epsIsDiagonal(const MxDimMatrix<double, 3> & epsTensor) const {
  double eps = std::numeric_limits<double>::epsilon();
  if (abs(epsTensor(0, 1)) > eps or abs(epsTensor(0, 2)) > eps or abs(epsTensor(1, 2)) > eps)
    return false;
  else
    return true;
}

template<size_t DIM>
MxDimMatrix<double, DIM> MxYeeMitInvEps<DIM>::aveInvEps(MxDimVector<double, DIM> pt, double vfrac,
    MxDimMatrix<double, DIM> epsIn, MxDimMatrix<double, DIM> epsOut) const {

  MxDimVector<double, DIM> normal = sim->getDielectric(0)->getShape()->normal(pt);
  MxDimMatrix<double, DIM> nn(normal, normal), eye(MxDimVector<double, DIM>(1));
  MxDimMatrix<double, DIM> nnC = eye - nn;
  MxDimMatrix<double, DIM> gammaIn, gammaOut, tauIn, tauOut, aveTau, aveEps;
  double tauNN;

  gammaIn = eye + nn * (eye - epsIn) / normal.dot(epsIn * normal);
  gammaOut = eye + nn * (eye - epsOut) / normal.dot(epsOut * normal);

  tauIn = gammaIn.transpose() * (nnC * epsIn * nnC - nn * epsIn * nn) *
    gammaIn;
  tauOut = gammaOut.transpose() * (nnC * epsOut * nnC - nn * epsOut * nn) *
    gammaOut;

  aveTau = vfrac * tauIn + (1.0 - vfrac) * tauOut;
  tauNN = normal.dot(aveTau * normal);

  aveEps = nnC * (aveTau - aveTau * nn * aveTau / tauNN) * nnC;
  aveEps -= nn / tauNN;
  aveEps -= nn * aveTau * nnC / tauNN;
  aveEps -= nnC * aveTau * nn / tauNN;
  return aveEps.inv();
}


template<>
MxPolytope<1> * MxYeeMitInvEps<1>::getVoxel() const {return NULL;}

template<>
MxPolytope<3> * MxYeeMitInvEps<3>::getVoxel() const {return new MxCartBox(grid->getCellSize());}

template<>
MxPolytope<2> * MxYeeMitInvEps<2>::getVoxel() const {
  return new MxCartRect<2>('z', grid->getCellSize()[0], grid->getCellSize()[1]);
}

template<size_t DIM>
void MxYeeMitInvEps<DIM>::fillAveEps() {
  // one guard cell domain:
  MxGridDomain<DIM> domain = grid->getGridDomain(1);

  MxPolytope<DIM> * voxel = getVoxel();

  MxDimMatrix<double, 3> epsIn = real(sim->getDielectric(0)->getEps());
  MxDimMatrix<double, 3> epsOut = real(sim->getBackgroundDielectric()->getEps());

  MxDimMatrix<double, DIM> epsDimIn, epsDimOut, invEpsIn, invEpsOut;
  for (size_t i = 0; i < DIM; ++i) {
    for (size_t j = 0; j < DIM; ++j) {
      epsDimIn(i, j) = epsIn(i, j);
      epsDimOut(i, j) = epsOut(i, j);
    }
  }
  invEpsIn = epsDimIn.inv();
  invEpsOut = epsDimOut.inv();

  size_t numFullCells = domain.getNumFullCells();

  epsOff.resize(numFullCells);
  for (int comp = 0; comp < DIM; ++comp)
    epsDiag.push_back(std::vector<MxDimMatrix<double, DIM> >(numFullCells));

  double vfrac;
  MxDimVector<int, DIM> cell;
  MxDimVector<double, DIM> loc;
  for (size_t i = 0; i < numFullCells; ++i) {
    cell = domain.fullIndxToCell(i);

    // get nodal average epsilon matrix
    loc = grid->nodeCoord(cell);
    vfrac = voxel->volumeFraction(*sim->getDielectric(0)->getShape(), loc);
    if (vfrac == 1.0)
      epsOff[i] = invEpsIn;
    else if (vfrac == 0.0)
      epsOff[i] = invEpsOut;
    else
      epsOff[i] = aveInvEps(loc, vfrac, epsDimIn, epsDimOut);

    // get average epsilons at component locations
    for (int comp = 0; comp < DIM; ++comp) {
      loc = efield->getCompCoord(comp, cell);
      vfrac = voxel->volumeFraction(*sim->getDielectric(0)->getShape(), loc);
      if (vfrac == 1.0)
        epsDiag[comp][i] = invEpsIn;
      else if (vfrac == 0.0)
        epsDiag[comp][i] = invEpsOut;
      else
        epsDiag[comp][i] = aveInvEps(loc, vfrac, epsDimIn, epsDimOut);
    }
  }

  delete voxel;
}

template<>
void MxYeeMitInvEps<3>::update(int comp0, MxDimVector<int, 3> cell,
    std::vector<double> & vals, std::vector<int> & inds, std::vector<double> &
    factors) {}

template<>
void MxYeeMitInvEps<1>::update(int comp0, MxDimVector<int, 1> cell,
    std::vector<double> & vals, std::vector<int> & inds, std::vector<double> &
    factors) {}

template<>
void MxYeeMitInvEps<2>::update(int comp0, MxDimVector<int, 2> cell,
    std::vector<double> & vals, std::vector<int> & inds, std::vector<double> &
    factors) {
  vals.clear();
  inds.clear();
  factors.clear();

  MxGridDomain<2> domain = grid->getGridDomain(1);

  int comp1 = (comp0 + 1) % 2;
  int fullIndx = domain.cellToFullIndx(cell);

  vals.push_back(epsDiag[comp0][fullIndx](comp0, comp0));
  inds.push_back(efield->globCompIndx(comp0, cell));
  if (this->imaginary)
    factors.push_back(efield->getCompFactor(comp0, cell).imag());
  else
    factors.push_back(efield->getCompFactor(comp0, cell).real());

  // 
  for (int i = 0; i < 2; ++i) {
    cell[comp0] += i;
    for (int j = 0; j < 2; ++j) {
      cell[comp1] -= i;
      fullIndx = domain.cellToFullIndx(cell);
      vals.push_back(0.25 * epsOff[fullIndx](comp0, comp1));
      inds.push_back(efield->globCompIndx(comp1, cell));
      if (this->imaginary)
        factors.push_back(efield->getCompFactor(comp1, cell).imag());
      else
        factors.push_back(efield->getCompFactor(comp1, cell).real());
      cell[comp1] += i;
    }
    cell[comp0] -= i;
  }
}


template<size_t DIM>
void MxYeeMitInvEps<DIM>::setMatrix2dTE3d() {
  MxDimMatrix<double, 3> epsIn = real(sim->getDielectric(0)->getEps());
  MxDimMatrix<double, 3> epsOut = real(sim->getBackgroundDielectric()->getEps());

  MxGridFieldIter<DIM> fieldIter(efield);

  int comp0;
  MxDimVector<int, DIM> cell;

  int row, col;
  double val;

  std::vector<double> factors, vals;
  std::vector<int> inds;

  // TE: Bz, Ex, Ey
  for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
    comp0 = fieldIter.getComp(); // Ex : 0; Ey : 1
    cell = fieldIter.getCell();
    row = fieldIter.getGlobCompIndx();

    update(comp0, cell, vals, inds, factors);

    for (size_t i = 0; i < vals.size(); ++i) {
      if (factors[i] != 0.0) {
        val = factors[i] * vals[i];
        col = inds[i];
        Epetra_CrsMatrix::InsertGlobalValues(row, 1, &val, &col);
      }
    }
  }

  Epetra_CrsMatrix::FillComplete(efield->getMap(), efield->getMap());
}

template<size_t DIM>
void MxYeeMitInvEps<DIM>::setMatrix2dTM() {
  MxDimMatrix<double, 3> epsIn = real(sim->getDielectric(0)->getEps());
  MxDimMatrix<double, 3> epsOut = real(sim->getBackgroundDielectric()->getEps());

  bool inIsDiag = epsIsDiagonal(epsIn);
  bool outIsDiag = epsIsDiagonal(epsOut);

  MxGridFieldIter<DIM> fieldIter(efield);

  MxDimVector<int, DIM> cell;

  int row;
  double val;
  double lfrac0, afrac0;

  // TM: Ez, Bx, By
  for (fieldIter.begin(); !fieldIter.atEnd(); fieldIter.bump()) {
    cell = fieldIter.getCell();
    row = efield->globCompIndx(0, cell);

    lfrac0 = efield->getCompFrac(0, cell, sim->getDielectric(0)->getName());
    afrac0 = dfield->getCompFrac(0, cell, sim->getDielectric(0)->getName()); 

    val = 1.0 / (afrac0 * epsIn(2, 2) + (1.0 - afrac0) * epsOut(2, 2));

    Epetra_CrsMatrix::InsertGlobalValues(row, 1, &val, &row);
  }

  Epetra_CrsMatrix::FillComplete(efield->getMap(), efield->getMap());
}

template<size_t DIM>
Epetra_CrsMatrix MxYeeMitInvEps<DIM>::cellAveInvEpsOperator() const {
  size_t numEComps = efield->getNumComps();

  MxDimMatrix<double, 3> epsIn = real(sim->getDielectric(0)->getEps());
  MxDimMatrix<double, 3> epsOut = real(sim->getBackgroundDielectric()->getEps());

  double invAveEpsIn, invAveEpsOut;
  if (DIM == 2 and numEComps == 1) {
    invAveEpsIn = 1.0 / epsIn(2, 2);
    invAveEpsOut = 1.0 / epsOut(2, 2);
  }
  else if (DIM == 2 and numEComps == 2) {
    invAveEpsIn = 2.0 / (epsIn(0, 0) + epsIn(1, 1));
    invAveEpsOut = 2.0 / (epsOut(0, 0) + epsOut(1, 1));
  }
  else {
    invAveEpsIn = 3.0 / epsIn.trace();
    invAveEpsOut = 3.0 / epsOut.trace();
  }

  const MxGridField<DIM> & psifield = sim->getField("psifield");
  Epetra_CrsMatrix res(Copy, psifield.getMap(), 1, true);

  MxGridFieldIter<DIM> iter(&psifield);
  double frac;
  int row;
  double val;

  for (iter.begin(); !iter.atEnd(); iter.bump()) {
    frac = psifield.getCompFrac(iter.getComp(), iter.getCell(), sim->getDielectric(0)->getName());
  
    row = iter.getGlobCompIndx();
    val = invAveEpsIn * frac + (1.0 - frac) * invAveEpsOut;

    res.InsertGlobalValues(row, 1, &val, &row);
  }

  res.FillComplete(psifield.getMap(), psifield.getMap());
  return res;
}


template class MxYeeMitInvEps<1>;
template class MxYeeMitInvEps<2>;
template class MxYeeMitInvEps<3>;
