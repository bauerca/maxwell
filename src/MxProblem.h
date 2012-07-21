#ifndef MX_PROBLEM
#define MX_PROBLEM

#include <string>
#include <set>
#include <vector>

#include "MxShapeDefs.h"
#include "MxGrid.h"
#include "MxDielectric.hpp"
#include "MxEMSim.h"
#include "MxPointCloud.h"
#include "MxSolver.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_Comm.h"

template<size_t DIM>
class MxProblem {
  public:
    MxProblem(const Teuchos::XMLObject & simXML, RCP<MxComm> theComm);

    ~MxProblem();

    void parseNode(const Teuchos::XMLObject & node);
    void preParse(const Teuchos::XMLObject & node);
    void postParse(const Teuchos::XMLObject & node);

    size_t getShape(const Teuchos::XMLObject & shapeNode);
    void getEigensolver(const Teuchos::XMLObject & node);
    void getLinearSolver(const Teuchos::XMLObject & node);
    void getGrid(const Teuchos::XMLObject & node);
    void getPEC(const Teuchos::XMLObject & node);
    void getPML(const Teuchos::XMLObject & node);
    void getDielectric(const Teuchos::XMLObject & node);
    void getMu(const Teuchos::XMLObject & node);
    void getBCs(const Teuchos::XMLObject & node);
    void getOutput(Teuchos::XMLObject const & node);

    void solve();

    static std::set<std::string> enumShapeTags();
    static const std::set<std::string> shapeTags;

    static int const numEigSolverParams;
    static std::string const eigSolverParams[][3];

    static int const numLinSolverParams;
    static std::string const linSolverParams[][3];

    // preconditioner stuff
    static int const numAMGParams;
    static std::string const amgParams[][3];

    static int const numILUParams;
    static std::string const iluParams[][3];


  private:
    
    void addParam(Teuchos::XMLObject const & node,
        std::string const & prefix, std::string const * set);

    int pid;

    std::string inFile;

    bool gridSet;

    bool eigSet;

    Teuchos::ParameterList pList;

    // container for ALL shapes
    std::vector<MxShape<DIM> *> shapes;

    std::vector<size_t> dielShapeIndices;
    std::vector<size_t> muShapeIndices;
    std::vector<size_t> pmlShapeIndices;
    size_t pecShapeIndx;

    Teuchos::RCP<MxGrid<DIM> > grid;
    Teuchos::RCP<MxDielectric<DIM> > diel;

    RCP<MxEMSim<DIM> > sim;

    //Teuchos::RCP<MxSolver<DIM> > solver;

    //std::vector<std::vector<MxDimVector<double, DIM> > > pointsSets;
    std::vector<MxPointCloud<DIM> > pointClouds;

    RCP<MxComm> comm;

};



#endif
