#ifndef MX_EM_SIM
#define MX_EM_SIM

#include <vector>
#include <string>

#include "MxGrid.h"
#include "MxGridField.hpp"
#include "MxDielectric.hpp"
#include "MxMu.hpp"
#include "MxPML.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

template<size_t DIM>
class MxEMSim {
  public:
    MxEMSim();

    MxEMSim(const MxGrid<DIM> * aGrid, Teuchos::ParameterList aPList);

    void setGrid(const MxGrid<DIM> * aGrid) {grid = aGrid;}

    const MxGrid<DIM> & getGrid() const {return *grid;}

    void setParameters(Teuchos::ParameterList aPList) {pList = aPList;}

    Teuchos::ParameterList getParameters() const {return pList;}

    /** set the background dielectric object.
     * This passes ownership to the EMSim object.
     */
    void setBackgroundDielectric(MxDielectric<DIM> const * dielectric) {
      mBGDiel = Teuchos::rcp(dielectric);
    }

    MxDielectric<DIM> const * getBackgroundDielectric() const {return mBGDiel.get();}

    /** add a dielectric to the EMSim's list of dielectrics.
     * This passes ownership to the EMSim object.
     */
    void addDielectric(MxDielectric<DIM> const * dielectric) {
      mDiels.push_back(Teuchos::rcp(dielectric));
      mNumDiels++;
      mHasDiel = true;
    }

    MxDielectric<DIM> const * getDielectric(size_t index) const {return mDiels[index].get();}

    bool hasDielectric() const {return mHasDiel;}

    size_t numDielectrics() const {return mNumDiels;}

    /** set the background mu object.
     * This passes ownership to the EMSim object.
     */
    void setBackgroundMu(MxMu<DIM> const * mu) {
      mBGMu = Teuchos::rcp(mu);
    }

    MxMu<DIM> const * getBackgroundMu() const {return mBGMu.get();}

    /** add a new mu material to the EMSim's list of mu material.
     * This passes ownership to the EMSim object.
     */
    void addMu(MxMu<DIM> const * mu) {
      mMus.push_back(Teuchos::rcp(mu));
      mNumMus++;
      mHasMu = true;
    }

    MxMu<DIM> const * getMu(size_t index) const {return mMus[index].get();}

    bool hasMu() const {return mHasMu;}

    size_t numMus() const {return mNumMus;}

    /** add a PML to the Sim's list of pmls. This passes ownership to the 
     * EMSim object.
     */
    void addPML(MxPML<DIM> const * pml) {
      mPmls.push_back(Teuchos::rcp(pml));
      mNumPmls++;
      mHasPml = true;
      mHasDiel = true;
      mHasMu = true;
    }

    bool hasPML() const {return mHasPml;}

    size_t numPMLs() const {return mNumPmls;}

    MxPML<DIM> const * getPML(size_t index) const {return mPmls[index].get();}

    void setPEC(const MxShape<DIM> * thePEC, double dmFrac) {
      pec = thePEC;
      mDMFrac = dmFrac;
    }

    const MxShape<DIM> * getPEC() const {return pec;}

    bool hasPEC() const {return pec != 0;}

    const MxGridField<DIM> & getField(std::string fieldName) const;

    int setup();

  private:
    const MxGrid<DIM> * grid;

    //// Dielectrics ////
    Teuchos::RCP<MxDielectric<DIM> const> mBGDiel;

    bool mHasDiel;

    int mNumDiels;

    // EMSim owns the dielectrics
    std::vector<Teuchos::RCP<MxDielectric<DIM> const> > mDiels;
    ////////////////////

    //// Mus ////
    Teuchos::RCP<MxMu<DIM> const> mBGMu;

    bool mHasMu;

    int mNumMus;

    // EMSim owns the mus
    std::vector<Teuchos::RCP<MxMu<DIM> const> > mMus;
    ////////////////////

    //// Pmls ////
    bool mHasPml;

    int mNumPmls;

    // EMSim owns the pmls
    std::vector<Teuchos::RCP<MxPML<DIM> const> > mPmls;
    //////////////

    const MxShape<DIM> * pec;

    double mDMFrac;

    Teuchos::ParameterList pList;

    std::vector<Teuchos::RCP<MxGridField<DIM> > > fields;

    std::vector<std::string> fieldNames;

};

#endif
