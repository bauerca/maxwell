
#include "MxEMSim.h"

#include "MxYeeFitEField.h"
#include "MxYeeFitDField.h"
#include "MxYeeFitBField.h"
#include "MxYeeFitHField.h"
#include "MxYeePsiField.h"

#include "MxTypes.h"

template<size_t DIM>
MxEMSim<DIM>::MxEMSim() : grid(0), pec(0),
mDMFrac(0),
mHasDiel(false), mNumDiels(0),
mHasMu(false), mNumMus(0),
mHasPml(false), mNumPmls(0),
mBGDiel(rcp(new MxDielectric<DIM>(
NULL, MxDimMatrix<MxComplex, 3>::I(), std::string("bg dielectric")))),
mBGMu(rcp(new MxMu<DIM>(
NULL, MxDimMatrix<MxComplex, 3>::I(), std::string("bg mu")))) {
  //setup();
}

template<size_t DIM>
MxEMSim<DIM>::MxEMSim(const MxGrid<DIM> * aGrid, Teuchos::ParameterList aPList) :
grid(aGrid), pList(aPList), pec(0),
mDMFrac(0),
mHasDiel(false), mNumDiels(0),
mHasMu(false), mNumMus(0),
mHasPml(false), mNumPmls(0),
mBGDiel(rcp(new MxDielectric<DIM>(
NULL, MxDimMatrix<MxComplex, 3>::I(), std::string("bg dielectric")))),
mBGMu(rcp(new MxMu<DIM>(
NULL, MxDimMatrix<MxComplex, 3>::I(), std::string("bg mu")))) {
  //setup();
}


template<size_t DIM>
const MxGridField<DIM> & MxEMSim<DIM>::getField(std::string fieldName) const {
  size_t numFields = fieldNames.size();
  for (size_t i = 0; i < numFields; ++i)
    if (fieldName == fieldNames[i])
      return *fields[i];

  std::cout << "MxEMSim<" << DIM << ">::getField: field '" << fieldName
            << "' not found.";
  throw 1;
}


template<size_t DIM>
int MxEMSim<DIM>::setup() {
  if (grid == 0) {
    std::cout << "MxEMSim<" << DIM << ">::setup(): grid not set.\n";
    throw 1;
  }

  int pid = grid->getComm()->myPID();

  fieldNames.clear();
  fields.clear();

  MxPolType pol = pList.get("polarization", TE);
  MxDimVector<MxBCType, DIM> lower =
      pList.get("lower bcs", MxDimVector<MxBCType, DIM>(PERIODIC));
  MxDimVector<MxBCType, DIM> upper =
      pList.get("upper bcs", MxDimVector<MxBCType, DIM>(PERIODIC));
  MxDimVector<double, DIM> phaseShifts = 
      pList.get("phase shifts", MxDimVector<double, DIM>(0.0));

  // setup fields
  // Yee Electric field
  fieldNames.push_back("efield");
  MxYeeFitEField<DIM> * efield = new MxYeeFitEField<DIM>(grid, pol);
  efield->setBCs(lower, upper);
  fields.push_back(rcp(efield));

  // Yee magnetic field
  fieldNames.push_back("bfield");
  MxYeeFitBField<DIM> * bfield = new MxYeeFitBField<DIM>(grid, pol);
  bfield->setDMFrac(mDMFrac);
  bfield->setBCs(lower, upper);
  fields.push_back(rcp(bfield));

  // Yee FIT electric displacement field
  if (this->hasDielectric() or this->hasPML()) {
    fieldNames.push_back("dfield");
    MxYeeFitDField<DIM> * dfield = new MxYeeFitDField<DIM>(grid, pol);
    dfield->setBCs(lower, upper);
    fields.push_back(rcp(dfield));
  }

  if (this->hasMu() or this->hasPML()) {
    fieldNames.push_back("hfield");
    MxYeeFitHField<DIM> * hfield = new MxYeeFitHField<DIM>(grid, pol);
    hfield->setBCs(lower, upper);
    fields.push_back(rcp(hfield));
  }

  // Yee FIT scalar field for divergence of B
  if (DIM == 3 or (DIM == 2 and pol == TM)) {
    fieldNames.push_back("psifield");
    MxYeePsiField<DIM> * psifield = new MxYeePsiField<DIM>(grid, bfield);
    psifield->setBCs(lower, upper);
    fields.push_back(rcp(psifield));
  }

  // set phase shifts on fields
  for (size_t i = 0; i < fields.size(); i++) {
    fields[i]->setPhaseShifts(phaseShifts);
  }

  // pass PECs to fields
  if (pec != 0) {
    if (pid == 0) std::cout << "Setting PEC on all fields...\n";
    for (size_t i = 0; i < fields.size(); i++) {
      if (pid == 0) std::cout << "  Setting PEC on " << fieldNames[i] << "...\n";
      // One guard cell. Set region.
      fields[i]->addShapeRep(*pec, "pec", 1, true);
    }
  }

  // pass dielectric objects to fields
  if (this->hasDielectric()) {
    if (pid == 0) std::cout << "Setting dielectrics on fields...\n";

    MxDielectric<DIM> const * diel;
    for (size_t i = 0; i < fields.size(); i++) {
      if (fieldNames[i] == "hfield" or fieldNames[i] == "bfield")
        continue;

      if (pid == 0) std::cout << "  Setting dielectrics on " << fieldNames[i] << "...\n";
      for (size_t j = 0; j < this->numDielectrics(); ++j) {
        diel = this->getDielectric(j);
        // One guard cell. Do not set region.
        fields[i]->addShapeRep(*diel->getShape(), diel->getName(), 1, false);
      }
    }
  }

  // pass mu objects to fields
  if (this->hasMu()) {
    if (pid == 0) std::cout << "Setting mu objects on fields...\n";

    MxMu<DIM> const * muObj;
    for (size_t i = 0; i < fields.size(); i++) {
      if (fieldNames[i] == "efield" or fieldNames[i] == "dfield")
        continue;

      if (pid == 0) std::cout << "  Setting mu objects on " << fieldNames[i] << "...\n";
      for (size_t j = 0; j < this->numMus(); ++j) {
        muObj = this->getMu(j);
        // One guard cell. Do not set region.
        fields[i]->addShapeRep(*muObj->getShape(), muObj->getName(), 1, false);
      }
    }
  }

// No need to set pml shapes on fields, since currently the FIT technique
// is not used on pmls
#if 0
  // pass PML objects to fields
  if (this->hasPML()) {
    if (pid == 0) std::cout << "Setting PMLs on all fields...\n";

    MxPML<DIM> const * pml;
    for (size_t i = 0; i < fields.size(); i++) {
      if (pid == 0) std::cout << "  Setting PMLs on " << fieldNames[i] << "...\n";
      for (size_t j = 0; j < this->numPMLs(); ++j) {
        pml = this->getPML(j);
        // One guard cell. Do not set region.
        fields[i]->addShapeRep(*pml->getShape(), pml->getName(), 1, false);
      }
    }
  }
#endif

  // if not pec, field maps must be set
  for (size_t i = 0; i < fields.size(); i++) {
    if (pec == 0 and (fieldNames[i] == "efield" or
        fieldNames[i] == "bfield" or
        fieldNames[i] == "psifield"))
      fields[i]->setMap();
  }
  // maps for D and H are always the same as E and B, resp.
  for (size_t i = 0; i < fields.size(); i++) {
    if (fieldNames[i] == "dfield")
      fields[i]->setMap(this->getField("efield"));
    else if (fieldNames[i] == "hfield")
      fields[i]->setMap(this->getField("bfield"));

    //std::cout << fieldNames[i] << " map addr: " << &(fields[i]->getMap()) << "\n";
    //if (pid == 0) std::cout << "  " << fieldNames[i] << ": " 
    //  << fields[i]->getMap().getGlobalNumElements() << " values\n";
  }

  
  //int nx = grid->getResolution()[0];
  //int ny = grid->getResolution()[1];
  //MxDimVector<int, DIM> cell; 
  //cell[0] = nx / 2;
  //cell[1] = ny / 2;
  //cell[2] = 2;
  //std::cout << "cell: "; cell.print();
  //std::cout << "Bx comp frac: " << fields[1]->getCompFrac(0, cell, "pec");
  //std::cout << "By comp frac: " << fields[1]->getCompFrac(1, cell, "pec");
  //std::cout << "Bz comp frac: " << fields[1]->getCompFrac(2, cell, "pec");
  //cell[2]++;
  //std::cout << "cell: "; cell.print();
  //std::cout << "Bx comp frac: " << fields[1]->getCompFrac(0, cell, "pec");
  //std::cout << "By comp frac: " << fields[1]->getCompFrac(1, cell, "pec");
  //std::cout << "Bz comp frac: " << fields[1]->getCompFrac(2, cell, "pec");
  //cell[2] -= 2;
  //std::cout << "cell: "; cell.print();
  //std::cout << "Bx comp frac: " << fields[1]->getCompFrac(0, cell, "pec");
  //std::cout << "By comp frac: " << fields[1]->getCompFrac(1, cell, "pec");
  //std::cout << "Bz comp frac: " << fields[1]->getCompFrac(2, cell, "pec");

  return 1;

}

template class MxEMSim<1>;
template class MxEMSim<2>;
template class MxEMSim<3>;
