
#include <cmath>

#include "MxGridField.hpp"
#include "MxGridFieldIter.hpp"
#include "MxUtil.hpp"

#include "hdf5.h"

#include "Epetra_Comm.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#endif

template<size_t DIM>
void MxGridField<DIM>::initBCs() {
  mCompLowerBCs.resize(numComps, MxDimVector<MxBCType, DIM>(PERIODIC));
  mCompUpperBCs.resize(numComps, MxDimVector<MxBCType, DIM>(PERIODIC));

  for (size_t i = 0; i < DIM; ++i)
    mPhaseFactors[i] = MxComplex(1.0, 0.0);
}

template<size_t DIM>
void MxGridField<DIM>::setCompBCs(size_t comp,
MxDimVector<MxBCType, DIM> lowerBCs,
MxDimVector<MxBCType, DIM> upperBCs) {
  mCompLowerBCs[comp] = lowerBCs;
  mCompUpperBCs[comp] = upperBCs;
}

template<size_t DIM>
void MxGridField<DIM>::setPhaseShifts(MxDimVector<double, DIM> phaseShifts) {
  MxComplex I(0.0, 1.0);
  for (size_t i = 0; i < DIM; ++i)
    mPhaseFactors[i] = exp(I * phaseShifts[i]);
}

template<size_t DIM>
bool MxGridField<DIM>::useCompInMap(size_t comp,
MxDimVector<int, DIM> cell) const {
  if (regionSet and
      this->getCompFrac(comp, cell, regionName) == 0.0)
    return false;

  MxDimVector<int, DIM> n = grid->getResolution();
  //std::cout << "ZERO: " << ZERO << ", PER: " << PERIODIC << "\n";

  MxBCType lbc, ubc;
  double xi;
  for (size_t i = 0; i < DIM; ++i) {
    lbc = mCompLowerBCs[comp][i];
    ubc = mCompUpperBCs[comp][i];
    xi = compCellCoords[comp][i];
    //std::cout << "celli " << cell[i] << ", xi " << xi << "\n";
    //std::cout << "celli == ni: " << (cell[i] == n[i]) << "\n";
    //std::cout << "xi == 0.0: " << (xi == 0.0) << "\n";
    //std::cout << "lbc: " << lbc << ", ubc: " << ubc << "\n";

    // case 1: on simulation lower bounds
    if (cell[i] == 0 and xi == 0.0) {
      if (lbc == ZERO)
        return false;
    }
    // case 2: above simulation upper bounds
    else if (cell[i] == n[i] and xi > 0.0) {
      return false;
    }
    // case 3: on simulation upper bounds
    else if (cell[i] == n[i] and xi == 0.0) {
      if (ubc == ZERO or ubc == PERIODIC)
        return false;
    }
  }
  return true;
}

template<size_t DIM>
MxDimVector<int, DIM> MxGridField<DIM>::getInteriorComp(size_t comp,
MxDimVector<int, DIM> cell) const {
  MxDimVector<int, DIM> newCell(cell);

  MxBCType lbc, ubc;
  double xi;

  MxDimVector<int, DIM> n = grid->getResolution();

  // for each dimension
  for (size_t i = 0; i < DIM; ++i) {
    lbc = mCompLowerBCs[comp][i];
    ubc = mCompUpperBCs[comp][i];
    xi = compCellCoords[comp][i];

    // case 1: component is below simulation lower bound and
    //         is on lower bound of grid cell
    if (cell[i] < 0 and xi == 0.0) {
      if (lbc == PERIODIC)
        newCell[i] = n[i] + cell[i];
      else
        newCell[i] = -cell[i];
    }
    // case 2: component is below simulation lower bound but
    //         somewhere inside the cell
    if (cell[i] < 0) {
      if (lbc == PERIODIC)
        newCell[i] = n[i] + cell[i];
      else
        newCell[i] = -cell[i] - 1;
    }
    // case 3: component is on simulation upper bound
    else if (cell[i] == n[i] and xi == 0.0) {
      if (ubc == PERIODIC)
        newCell[i] = 0;
    }
    // case 4: component is above simulation upper bound
    //         and is on cell lower bound
    else if (cell[i] >= n[i] and xi == 0.0) {
      if (ubc == PERIODIC)
        newCell[i] = cell[i] - n[i];
      else
        newCell[i] = n[i] - (cell[i] - n[i]);
    }
    // case 5: component is above simulation upper bound
    //         and is somewhere inside cell
    else if (cell[i] >= n[i]) {
      if (ubc == PERIODIC)
        newCell[i] = cell[i] - n[i];
      else
        newCell[i] = n[i] - (cell[i] - n[i] + 1);
    }

  }

  //if (not (newCell == cell)) {
  //  std::cout << "Field " << fieldName << ", comp = " << comp << "\n";
  //  std::cout << "  Interiorized cell: " << newCell;
  //  std::cout << "  Original cell: " << cell;
  //}

  return newCell;
}

template<size_t DIM>
MxComplex MxGridField<DIM>::getCompFactor(size_t comp,
MxDimVector<int, DIM> cell) const {
  if (regionSet)
    if (this->getCompFrac(comp, cell, regionName) == 0.0)
      return 0.0;

  MxComplex res(1.0, 0.0);

  MxBCType lbc, ubc;
  double xi;
  MxDimVector<int, DIM> n = grid->getResolution();

  for (size_t i = 0; i < DIM; ++i) {
    lbc = mCompLowerBCs[comp][i];
    ubc = mCompUpperBCs[comp][i];
    xi = compCellCoords[comp][i];

    // case 1: component is below simulation lower bound
    if (cell[i] < 0) {
      if (lbc == PERIODIC)
        res /= mPhaseFactors[i];
      else if (lbc == ZERO)
        res *= -1.0;
    }
    // case 2: component is on simulation lower bound
    else if (cell[i] == 0 and xi == 0.0) {
      if (lbc == ZERO)
        res *= 0.0;
    }
    // case 3: component is on simulation upper bound
    else if (cell[i] == n[i] and xi == 0.0) {
      if (ubc == PERIODIC)
        res *= mPhaseFactors[i];
      else if (ubc == ZERO)
        res *= 0.0;
    }
    // case 4: component is above simulation upper bound
    else if (cell[i] >= n[i]) {
      if (ubc == PERIODIC)
        res *= mPhaseFactors[i];
      else if (ubc == ZERO)
        res *= -1.0;
    }
  }
  return res;
}

template<size_t DIM>
void MxGridField<DIM>::addShapeRep(const MxShape<DIM> & aShape, std::string shapeName, int numGuard, bool setRegion) {

  //shapeRepsDomains[name] = grid->getGridDomain(numGuard);
  shapeRepsDomains.insert(std::make_pair(shapeName, grid->getGridDomain(numGuard)));

  const MxGridDomain<DIM> & domain = shapeRepsDomains.find(shapeName)->second;

  size_t numFullCells = domain.getNumFullCells();

  std::vector<double> shapeRep(numComps * numFullCells);

  size_t fullCompIndx = 0;
  MxDimVector<double, DIM> p, node;
  MxDimVector<int, DIM> cell;
  for (size_t i = 0; i < numFullCells; i++) {
    cell = domain.fullIndxToCell(i);
    //node = grid->nodeCoord(cell);
    for (size_t comp = 0; comp < numComps; comp++) {
      shapeRep[fullCompIndx] = this->calcCompFrac(comp, cell, aShape);
      //p = node + compCellCoords[comp];
      //std::cout << "comp " << comp << " cell "; cell.print();
      //shapeRep[fullCompIndx] = compPtopes[comp]->volumeFraction(aShape, p);
      //std::cout << "  " << shapeRep[fullCompIndx] << "\n";
      fullCompIndx++;
    }
  }

  shapeReps[shapeName] = shapeRep;

  if (setRegion)
    MxGridField<DIM>::setRegion(shapeName);

  //MxUtil::printStdVector(compVolFracs);
}

template<size_t DIM>
void MxGridField<DIM>::setRegion(std::string shapeName) {
  if (!hasShapeRep(shapeName)) {
    std::cout << "MxGridField<" << DIM << ">::setRegion: shapeRep with name '" 
              << shapeName << "' does not exist.\n";
    throw 1;
  }
  else {
    regionName = shapeName;
    regionSet = true;
    this->setMap();
  }
}

template<size_t DIM>
void MxGridField<DIM>::setMap(const MxGridField<DIM> & partnerField) {
  if (partnerField.numComps != numComps) {
    std::cout << "MxGridField<DIM>::setMap(field): fields "
              << fieldName << " and " << partnerField.fieldName 
              << " do not have the same number of components, can't cross-set map.\n";
    throw 1;
  }
  
  mMap = partnerField.mMap;
  mapSet = true;
}

template<size_t DIM>
void MxGridField<DIM>::setMap() {

  MxGridDomain<DIM> * domain;
  if (regionSet)
    domain = &(shapeRepsDomains.find(regionName)->second);
  else
    domain = new MxGridDomain<DIM>(grid->getGridDomain(0));

#if 0
  // using a grid domain iterator
  MxGridDomainIter<DIM> gdIter = domain->interiorBegin();
  MxGridDomainIter<DIM> gdIterEnd = domain->interiorEnd();
  for (gdIter; gdIter != gdIterEnd; gdIter.interiorBump()) {
    cell = gdIter.getCell();
    fullCellIndx = gdIter.getFullDomainIndx();
  }
#endif

  int numIntCells = domain->getNumInteriorCells();
  int numIntComps = numComps * numIntCells;

  std::vector<MxIndex> globCompInds;
  globCompInds.reserve(numIntComps);

  MxDimVector<int, DIM> cell;
  MxIndex globCellIndx, globCompIndx;
  for (int i = 0; i < numIntCells; i++) {
    cell = domain->interiorIndxToCell(i);
    globCellIndx = grid->cellToGlobalIndx(cell);
    for (size_t comp = 0; comp < numComps; ++comp) {
      if (this->useCompInMap(comp, cell)) {
        globCompIndx = comp + numComps * globCellIndx;
        globCompInds.push_back(globCompIndx);
      }
    }
  }

  mMap = rcp(new MxMap(globCompInds, grid->getComm()));
  mapSet = true;

  if (!regionSet) delete domain;
}

template<size_t DIM>
bool MxGridField<DIM>::compInMap(size_t comp, MxDimVector<int, DIM> cell) const {
  MxIndex globIndx = this->globCompIndx(comp, cell);

  // first check if on local processor
  int pid;
  MxIndex lid;
  int res;
  //if (map->LID(globIndx) != -1)
  if (mMap->getLocalIndex(globIndx) != MxInvalidIndex)
    return true;
  else {
    // have to check globally
    //mMap->RemoteIDList(1, &globIndx, &pid, &lid);
    res = mMap->getRemoteIndexList(1, &globIndx, &pid, &lid);
    //if (lid != -1)
    if (res == -1)
      return false;
    else
      return true;
  }

}

#if 0
// old way. Faster?
template<size_t DIM>
void MxGridField<DIM>::setEpetraMap() {

  std::vector<int> globCompInds;

  MxDimVector<int, DIM> cell;
  size_t globCellIndx;
  int numIntCells;
  MxDimVector<int, DIM> gridRes = this->grid->getResolution();
  
  //MxGridDomainIter<DIM> iter;
  // if there is a restricted region set, then the MxMap should only
  // contain those field components that are partially or fully inside the region.
  if (regionSet) {
    const MxGridDomain<DIM> & domain = shapeRepsDomains.find(regionName)->second;
    //MxGridDomain<DIM> domain(grid->getGridDomain(shapeRepsGuardCells[region]));
    numIntCells = domain.getNumInteriorCells();

    const std::vector<double> & shapeRep = shapeReps[regionName];

    globCompInds.reserve(numComps * numIntCells);
    size_t fullCompIndx, fullCellIndx;
    for (int i = 0; i < numIntCells; i++) {
      cell = domain.interiorIndxToCell(i);
      // by default, do not include extra layer of cells at upper boundaries
      if (anyEqual(cell, gridRes))
        continue;
      fullCellIndx = domain.cellToFullIndx(cell);
      globCellIndx = grid->cellToGlobalIndx(cell);
// do this if you want to exclude all components lying
// outside of the region
#if 1
      for (size_t comp = 0; comp < numComps; comp++) {
        fullCompIndx = comp + numComps * fullCellIndx;
        if (shapeRep[fullCompIndx] > 0)
          globCompInds.push_back(comp + numComps * globCellIndx);
      }
#endif
// do this if you want to exclude only those components belonging
// to a cell that is completely outside of the region
#if 0
      bool anyInside = false;
      for (size_t comp = 0; comp < numComps; comp++) {
        fullCompIndx = comp + numComps * fullCellIndx;
        if (shapeRep[fullCompIndx] > 0) anyInside = true;
      }
      if (anyInside)
        for (size_t comp = 0; comp < numComps; comp++)
          globCompInds.push_back(comp + numComps * globCellIndx);
#endif
    }
  }
  else {
    // get this processor's grid domain w/ 0 guard cells
    MxGridDomain<DIM> domain(grid->getGridDomain(0));
    numIntCells = domain.getNumInteriorCells();

    globCompInds.resize(numComps * numIntCells);
    size_t intCompIndx = 0;
    for (int i = 0; i < numIntCells; i++) {
      cell = domain.interiorIndxToCell(i);
      globCellIndx = grid->cellToGlobalIndx(cell);
      for (size_t comp = 0; comp < numComps; comp++) {
        globCompInds[intCompIndx] = comp + numComps * globCellIndx;
        intCompIndx++;
      }
    }
  }

  map = Teuchos::rcp(new MxMap(-1, globCompInds.size(), &globCompInds[0], 0, grid->getComm()));

  mapSet = true;
}
#endif

template<size_t DIM>
RCP<MxVector<double> > MxGridField<DIM>::getAllCompFracs(
std::string shapeName) const {
  if (mMap == Teuchos::null) {
    std::cout << "MxGridField<" << DIM << ">::getAllCompFracs: need to call setMap() first!\n";
    throw 1;
  }
  if (!hasShapeRep(shapeName)) {
    std::cout << "MxGridField<" << DIM << ">::getAllCompFracs: shapeRep, '" << shapeName << "', does not exist.\n";
    throw 1;
  }

  //Epetra_Vector res(*map);
  RCP<MxVector<double> > res = rcp(new MxVector<double>(mMap));

  const MxGridDomain<DIM> & domain = shapeRepsDomains.find(shapeName)->second;
  //MxGridDomain<DIM> domain(grid->getGridDomain(shapeRepsGuardCells[region]));
  size_t numIntCells = domain.getNumInteriorCells();

  MxDimVector<int, DIM> cell;
  MxIndex localMapIndx = 0;
  MxIndex fullCellIndx;

  double frac;
  for (size_t i = 0; i < numIntCells; i++) {
    cell = domain.interiorIndxToCell(i);
    fullCellIndx = domain.cellToFullIndx(cell);
    for (size_t comp = 0; comp < numComps; comp++) {
      frac = shapeReps.find(shapeName)->second[comp + numComps * fullCellIndx];
      if (this->useCompInMap(comp, cell)) {
        res->replaceLocalValue(localMapIndx, frac);
        localMapIndx++;
      }
    }
  }

  return res;
}


template<size_t DIM>
template<typename Scalar>
RCP<MxMultiVector<Scalar> > MxGridField<DIM>::uniformFields(
MxGridField<DIM> const & field) {

  RCP<MxMultiVector<Scalar> > res;
  if (ScalarTraits<Scalar>::isComplex)
    res = rcp(new MxMultiVector<Scalar>(field.getMap(), 2 * field.getNumComps()));
  else
    res = rcp(new MxMultiVector<Scalar>(field.getMap(), field.getNumComps()));

  MxGridFieldIter<DIM> iter(&field);

  size_t comp;
  MxIndex row;
  MxDimVector<int, DIM> cell;
  Scalar val;
  for (iter.begin(); !iter.atEnd(); iter.bump()) {
    comp = iter.getComp();
    row = iter.getGlobCompIndx();
    cell = iter.getCell();
    if (field.regionSet) {
      if (field.getCompFrac(comp, cell, field.regionName) == 0.0) {
        continue;
      }
    }

    if (ScalarTraits<Scalar>::isComplex) {
      res->replaceGlobalValue(row, 2 * comp, ScalarTraits<Scalar>::one());
      res->replaceGlobalValue(row, 2 * comp + 1, MxUtil::i<Scalar>());
    }
    else
      res->replaceGlobalValue(row, comp, ScalarTraits<Scalar>::one());
  }

  res->normalize();
  return res;
}

template RCP<MxMultiVector<double> > MxGridField<1>::uniformFields<double>(
MxGridField<1> const & field);
template RCP<MxMultiVector<double> > MxGridField<2>::uniformFields<double>(
MxGridField<2> const & field);
template RCP<MxMultiVector<double> > MxGridField<3>::uniformFields<double>(
MxGridField<3> const & field);
template RCP<MxMultiVector<MxComplex> > MxGridField<1>::uniformFields<MxComplex>(
MxGridField<1> const & field);
template RCP<MxMultiVector<MxComplex> > MxGridField<2>::uniformFields<MxComplex>(
MxGridField<2> const & field);
template RCP<MxMultiVector<MxComplex> > MxGridField<3>::uniformFields<MxComplex>(
MxGridField<3> const & field);

template<size_t DIM>
template<typename Scalar>
void MxGridField<DIM>::maxValueLocation(MxMultiVector<Scalar> const & data,
int vec,
MxDimVector<int, DIM> & outCell, int & outComp, Scalar & outValue) const {
  
  MxGridFieldIter<DIM> iter(this);

  outValue = Teuchos::ScalarTraits<Scalar>::zero();
  Scalar val = Teuchos::ScalarTraits<Scalar>::zero();

  MxIndex row;
  MxDimVector<int, DIM> cell;
  double frac0 = 1;
  double frac1 = 1;
  int count = 0;
  for (iter.begin(); !iter.atEnd(); iter.bump()) {
    val = data(iter.getLocalCompIndx(), vec);
    if (abs(val) > abs(outValue)) {
      outValue = val;
      outCell = iter.getCell();
      outComp = iter.getComp();
      if (this->regionSet) {
        frac0 = this->getCompFrac(outComp, outCell, this->regionName);
        frac1 = this->getCompFrac((outComp + 1) % this->getNumComps(), outCell, this->regionName);
      }
    }
  }

  std::cout << "val=" << outValue << ", comp=" << outComp << ", cell="
      << outCell;
  std::cout << "  comp="<<outComp<<", frac=" << frac0 << "\n";
  std::cout << "  comp="<<((outComp+1)%this->getNumComps())<<", frac=" << frac1 << "\n";
}

template void MxGridField<1>::maxValueLocation<double>(
MxMultiVector<double> const & data, int vec,
MxDimVector<int, 1> & outCell, int & outComp, double & outValue) const;
template void MxGridField<2>::maxValueLocation<double>(
MxMultiVector<double> const & data, int vec,
MxDimVector<int, 2> & outCell, int & outComp, double & outValue) const;
template void MxGridField<3>::maxValueLocation<double>(
MxMultiVector<double> const & data, int vec,
MxDimVector<int, 3> & outCell, int & outComp, double & outValue) const;
template void MxGridField<1>::maxValueLocation<MxComplex>(
MxMultiVector<MxComplex> const & data, int vec,
MxDimVector<int, 1> & outCell, int & outComp, MxComplex & outValue) const;
template void MxGridField<2>::maxValueLocation<MxComplex>(
MxMultiVector<MxComplex> const & data, int vec,
MxDimVector<int, 2> & outCell, int & outComp, MxComplex & outValue) const;
template void MxGridField<3>::maxValueLocation<MxComplex>(
MxMultiVector<MxComplex> const & data, int vec,
MxDimVector<int, 3> & outCell, int & outComp, MxComplex & outValue) const;

template<size_t DIM>
template<typename Scalar>
void MxGridField<DIM>::zeroUnusedComponents(MxMultiVector<Scalar> & mv,
MxGridField<DIM> const & field) {
  if (!field.regionSet) {
    return;
  }

  MxGridFieldIter<DIM> iter(&field);

  size_t comp;
  MxIndex row;
  Scalar val;
  MxDimVector<int, DIM> cell;
  int count = 0;
  for (iter.begin(); !iter.atEnd(); iter.bump()) {
    comp = iter.getComp();
    cell = iter.getCell();
    row = iter.getGlobCompIndx();

    if (field.getCompFrac(comp, cell, field.regionName) == 0.0) {
      count++;
      for (int i = 0; i < mv.getNumVecs(); ++i) {
        //std::cout << "zeroing\n";
        mv.replaceGlobalValue(row, i, ScalarTraits<Scalar>::zero());
      }
    }
  }

  std::cout << count << " components zeroed\n";
}

template void MxGridField<1>::zeroUnusedComponents<double>(
MxMultiVector<double> & mv, MxGridField<1> const & field);
template void MxGridField<2>::zeroUnusedComponents<double>(
MxMultiVector<double> & mv, MxGridField<2> const & field);
template void MxGridField<3>::zeroUnusedComponents<double>(
MxMultiVector<double> & mv, MxGridField<3> const & field);
template void MxGridField<1>::zeroUnusedComponents<MxComplex>(
MxMultiVector<MxComplex> & mv, MxGridField<1> const & field);
template void MxGridField<2>::zeroUnusedComponents<MxComplex>(
MxMultiVector<MxComplex> & mv, MxGridField<2> const & field);
template void MxGridField<3>::zeroUnusedComponents<MxComplex>(
MxMultiVector<MxComplex> & mv, MxGridField<3> const & field);

template<size_t DIM>
RCP<MxMultiVector<double> > MxGridField<DIM>::coords() const {
  RCP<MxMultiVector<double> > res;
  res = rcp(new MxMultiVector<double>(mMap, DIM));

  MxGridFieldIter<DIM> iter(this);

  MxIndex row;
  MxDimVector<double, DIM> coord;
  for (iter.begin(); !iter.atEnd(); iter.bump()) {
    row = iter.getGlobCompIndx();
    coord = iter.getCoord();
    for (size_t i = 0; i < DIM; ++i)
      res->replaceGlobalValue(row, i, coord[i]);
  }

  return res;
}


#if 0
template<size_t DIM>
MxMap MxGridField<DIM>::createPointsMap(std::vector<MxDimVector<double, DIM> > const & points) const {
  MxGridDomain<DIM> domain = grid->getGridDomain(0);
  MxDimVector<double, DIM> lb = grid->nodeCoord(domain.getLowerBoundCell());
  MxDimVector<double, DIM> ub = grid->nodeCoord(domain.getUpperBoundCell());

  std::vector<MxDimVector<double, DIM> >::const_iterator iter;
  bool ptInDomain;
  for (iter = points.begin(); iter != points.end(); ++iter) {
    ptInDomain = true;
    for (size_t i = 0; i < DIM; ++i) {
      if ((*iter)[i] < lb[i] or (*iter)[i] > ub[i]) {
        ptInDomain = false;
        break;
      }
    }
    

  }



}
#endif



#if 0
template<size_t DIM>
void MxGridField<DIM>::load(std::string filename) {
  MxIO<DIM> io(&grid->getComm());

  MxHDF5<DIM> h5(&grid->getComm());

  h5.openFile(filename, 'r');

  char strbuf1[200], strbuf2[200];
  double * fracs;

  //if (h5.hasDSet("/pec/comp_fracs")) {
  if (h5.numDSets("/pec") > 0) {
    MxGridFieldData<DIM> * fracs;
    h5.getDomainDSet("pec/pec_comp_fracs", fracs);
    shapeReps["pec"] = Teuchos::rcp(fracs);
    shapeReps["pec"] = Teuchos::rcp(new Epetra_Vector(Copy, grid->fieldMap(this->numComps), fracs));
    delete[] fracs;
  }

  if (h5.numDSets("/dielectric") > 0) {
    for (int i = 0; i < h5.numDSets("/dielectric"); ++i) {
      sprintf(strbuf1, "/dielectric/diel%i_comp_fracs", i);
      h5.getDomainDSet(strbuf1, fracs);
      strbuf2 = h5.getDSetName(strbuf1).c_str();
      shapeReps[strbuf2] = Teuchos::rcp(new Epetra_Vector(Copy, grid->fieldMap(this->numComps), fracs));
      delete[] fracs;
    }
  }

}

template<size_t DIM>
void MxGridField<DIM>::save(std::string filename) const {
  MxHDF5<DIM> h5(&grid->getComm());
  h5.openFile(filename, 'w');

  typename std::map<std::string, Teuchos::RCP<MxGridFieldData<DIM> > >::const_iterator iter;
  for (iter = shapeReps.begin(); iter != shapeReps.end(); ++iter) {
    sprintf(buf, "/comp_fracs/%s", iter->first.c_str());
    h5.write(buf, *iter->second);
    h5.write(path, name, 
  }



}

#endif



template class MxGridField<1>;
template class MxGridField<2>;
template class MxGridField<3>;
