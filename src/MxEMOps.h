#ifndef MX_EM_OPS
#define MX_EM_OPS

#include <vector>
#include <string>

#include "MxTypes.h"
#include "MxGrid.h"
#include "MxGridField.hpp"
#include "MxCrsMatrix.hpp"
#include "MxEMSim.h"

#include "Teuchos_ParameterList.hpp"

template<size_t DIM, typename Scalar>
class MxEMOps {
  public:
    MxEMOps();

    //~MxEMOps();

    void printRCPs();

    MxEMOps(RCP<MxEMSim<DIM> > sim, bool initSingleOps = true);

    void setSim(RCP<MxEMSim<DIM> > sim) {mSim = sim;}

    /** 
     * Constructs single operators automatically based on simulation
     * requirements. Operators constructed are (if necessary):
     *
     *    curlB, curlE, divB, gradPsi, invEps, mu, dmL, dmA, dmVInv,
     *    invMuCellAve, invEpsCellAve
     *
     */
    void setSingleOps(RCP<MxEMSim<DIM> > sim);

    void setSingleOps() {setSingleOps(mSim);}

    /**
     * name must be one of:
     *
     *  curlE, curlB, divB, gradPsi, invEps, mu
     *
     */
    void setSingleOp(RCP<MxEMSim<DIM> > sim, std::string name) {};

    /**
     * Less typing if MxEMSim is already set for this object. Invokes
     * the above with the previously-provided MxEMSim. Will crash if
     * MxEMSim hasn't already been given.
     *
     */
    void setSingleOp(std::string name) {setSingleOp(mSim, name);}

    RCP<MxCrsMatrix<Scalar> > getOp(std::string name) const;
    
    void gatherOps(std::vector<std::string> const & names,
      std::vector<RCP<MxCrsMatrix<Scalar> > > & mats) const;

    /**
     * Make a new operator by multiplying existing operators. 
     * Domain of the first operator in the list is the domain of the
     * resultant operator; range of result is range of last operator in list.
     * Optionally stores result under 'newName.'
     *
     * Several things could fail: operator in list could not exist,
     * operators
     * may not have compatible domain/range maps.
     *
     */
    RCP<MxCrsMatrix<Scalar> > multiplyOps(std::vector<std::string> names,
      std::vector<bool> transposes, bool fillComplete, bool store = true,
      std::string name = std::string(""));

    /**
     * Make a new operator by adding existing operators. 
     * Domain and ranges of all operators in list must be identical.
     * Stores result under 'newName.'
     *
     * Several things could fail: operator in list could not exist,
     * operators
     * may not have compatible domain/range maps.
     *
     */
    RCP<MxCrsMatrix<Scalar> > addOps(std::vector<std::string> names,
      std::vector<Scalar> coeffs, bool store, std::string name);

    /**
    Saves all operators to MatrixMarket files.
    */
    void saveOps() const;

    void setMagVecLapl();

    void setPsiScaLapl();

    Teuchos::RCP<Epetra_CrsMatrix> getMagVecLapl() const;

  private:
    int pid;

    const MxGrid<DIM> * grid;

    Teuchos::ParameterList pList;

    RCP<MxEMSim<DIM> > mSim;

    typename std::map<std::string, RCP<MxCrsMatrix<Scalar> > > mOps;

};

#endif
