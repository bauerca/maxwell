#ifndef MX_GEO_MG_PREC
#define MX_GEO_MG_PREC

#include <vector>

#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_Preconditioner.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

class Epetra_Comm;
class MxMap;
class Epetra_MultiVector;

class MxGeoMultigridPrec : virtual public Epetra_Operator {

  public:
    MxGeoMultigridPrec() {}

    ~MxGeoMultigridPrec();

    MxGeoMultigridPrec(Teuchos::ParameterList mgPList);

    void setOperators(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs);

    void setCoarseners(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs);

    void setLHSCoarseners(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs);

    void setRHSCoarseners(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs);

    void setRefiners(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs);

    void setLHSRefiners(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs);

    void setRHSRefiners(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs);

    void setup();

    int ApplyInverse(const Epetra_MultiVector & b, Epetra_MultiVector & x) const;

    int Apply(const Epetra_MultiVector & x, Epetra_MultiVector & b) const {return 0;}

    const char * Label () const {
      return "Geometric multigrid preconditioner";
    }

    int SetUseTranspose (bool sut) {return 0;}

    bool UseTranspose () const {return false;}
    
    bool HasNormInf() const {return false;}

    double NormInf() const {return 0;}

    const Epetra_Comm & Comm() const {return ops[0]->Comm();}

    const MxMap & OperatorDomainMap() const {return ops[0]->DomainMap();}

    const MxMap & OperatorRangeMap() const {return ops[0]->RangeMap();}

  private:

    Teuchos::ParameterList pList;

    bool opsSet;

    bool cSet;

    bool rSet;

    bool rcf;
    
    //std::vector<const Epetra_CrsMatrix *> ops, coarseners, refiners;
    std::vector<Teuchos::RCP<Epetra_CrsMatrix> > ops, lhsCoarseners, rhsCoarseners, lhsRefiners, rhsRefiners;

    std::vector<Teuchos::RCP<Ifpack_Preconditioner> > smoothers;

    Teuchos::RCP<Epetra_CrsMatrix> fineOpCopy;

    std::vector<double> maxEigs;

    int levels;

    int cycles;
    
    int sweeps;
    
    int output;

    void removeConstField(Epetra_MultiVector & x) const;

    void smoothInterpolation(int level, Epetra_MultiVector & x) const;

    double getResidual(const Epetra_CrsMatrix & op, const Epetra_MultiVector & x, const Epetra_MultiVector & b) const;

    double constPart(const Epetra_MultiVector & x) const;

#if 0
    void GetBoundaryBulkResiduals(const Epetra_CrsMatrix & op, 
                                    const Epetra_MultiVector & x, 
                                    const Epetra_MultiVector & b,
                                    double & bndryRes,
                                    double & bulkRes,
                                    double & totalRes,
                                    int level) const;
#endif
    void refine(int targetLevel, const Epetra_MultiVector & coarse, Epetra_MultiVector & fine) const;

    void coarsen(int targetLevel, const Epetra_MultiVector & fine, Epetra_MultiVector & coarse) const;

    int vCycle(std::vector<Epetra_MultiVector*> & bvecs, std::vector<Epetra_MultiVector*> & xvecs, int startlevel, int cycle) const;

    int fullVCycle(std::vector<Epetra_MultiVector*> & bvecs, std::vector<Epetra_MultiVector*> & xvecs) const;

    void spaces(int level) const;

    mutable double totalTime;

    mutable int numApplyInverse;
};

#endif
