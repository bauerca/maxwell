
#include "MxProblem.h"
#include "MxShapeDefs.h"
#include "MxTypes.h"
#include "MxIO.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FileInputSource.hpp"

template<size_t DIM>
std::set<std::string> MxProblem<DIM>::enumShapeTags() {
  std::string data[] = {"Sphere", "Ellipsoid", "Torus", "Cylinder",
                        "Slab", "HalfSpace", "Cone", "Intersection",
                        "Difference", "Union", "Subtract"};
  return std::set<std::string>(data, data + 11);
}

template<size_t DIM>
const std::set<std::string> MxProblem<DIM>::shapeTags = MxProblem<DIM>::enumShapeTags();

template<size_t DIM>
int const MxProblem<DIM>::numEigSolverParams = 10;

// format: name, default, type
template<size_t DIM>
std::string const MxProblem<DIM>::eigSolverParams[numEigSolverParams][3] = {
    {"type", "krylov-schur", "string"},
    {"basis", "20", "int"},
    {"max restarts", "100", "int"},
    {"invert", "true", "bool"},
    {"block size", "1", "int"},
    {"nev", "10", "int"},
    {"tol", "1.0e-6", "double"},
    {"spectrum", "LM", "string"},
    {"output", "2", "int"},
    {"shift", "0.0", "double"}};


template<size_t DIM>
int const MxProblem<DIM>::numLinSolverParams = 3;

template<size_t DIM>
std::string const MxProblem<DIM>::linSolverParams[numLinSolverParams][3] = {
    {"type", "gmres", "string"},
    {"tol", "1.0e-6", "double"},
    {"basis", "40", "int"}};



template<size_t DIM>
int const MxProblem<DIM>::numAMGParams = 4;

template<size_t DIM>
std::string const MxProblem<DIM>::amgParams[numAMGParams][3] = {
    {"smoother", "Chebyshev", "string"},
    {"sweeps", "2", "int"},
    {"coarse smoother", "Chebyshev", "string"},
    {"max levels", "10", "int"}};


template<size_t DIM>
int const MxProblem<DIM>::numILUParams = 2;

template<size_t DIM>
std::string const MxProblem<DIM>::iluParams[numILUParams][3] = {
    {"drop tol", "1.e-3", "double"},
    {"fill", "10.", "double"}};



//// format: name, default, type
//template<size_t DIM>
//std::string const MxProblem<DIM>::linSolverParams[numLinSolverParams][3] = {
//    {"type", "krylov-schur", "string"},
//    {"basis", "20", "int"},
//    {"max restarts", "100", "string"},
//    {"invert", "true", "bool"},
//    {"block size", "1", "int"},
//    {"nev", "10", "int"},
//    {"tol", "1.0e-6", "double"},
//    {"spectrum", "LM", "string"},
//    {"output", "2", "int"},
//    {"shift", "0.0", "double"};

template<size_t DIM>
void MxProblem<DIM>::addParam(Teuchos::XMLObject const & node,
std::string const & prefix, std::string const * set) {
  std::string name = set[0];
  std::string dflt = set[1];
  std::string type = set[2];
  std::string sval = MxUtil::XML::getAttr(name, node, dflt);

  if (type == "int") {
    int val = atoi(sval.c_str());
    pList.set(prefix + name, val);
  }
  else if (type == "double") {
    double val = atof(sval.c_str());
    pList.set(prefix + name, val);
  }
  else if (type == "bool") {
    bool val = (sval == "true" ? true : false);
    pList.set(prefix + name, val);
  }
  else {
    pList.set(prefix + name, sval);
  }
}

template<size_t DIM>
MxProblem<DIM>::MxProblem(const Teuchos::XMLObject & simXML, RCP<MxComm> theComm) : 
sim(rcp(new MxEMSim<DIM>())), gridSet(false), eigSet(false), comm(theComm),
pid(theComm->myPID()) {
  std::cout << "Running simulation: " << simXML.getTag() << "\n";

  // set polarization
  std::string pol = MxUtil::XML::getAttr("polarization", simXML);
  if (pol == "TE")
    pList.set("polarization", TE);
  else if (pol == "TM")
    pList.set("polarization", TM);
  else {
    if (pid == 0) std::cout << "Invalid polarization: '" << pol << "'. Should be 'TE' or 'TM'.\n";
    exit(EXIT_FAILURE);
  }
  pList.set("is complex", false);
  pList.set("has curl null", DIM == 3 or (DIM == 2 and pol == "TM"));

  for (int i = 0; i < simXML.numChildren(); ++i)
    //parseNode(simXML.getChild(i));
    preParse(simXML.getChild(i));

  if (!eigSet)
    getEigensolver(Teuchos::XMLObject()); // set eigensolver with default params

  if (!gridSet)
    getGrid(Teuchos::XMLObject()); // set grid with default params

  sim->setGrid(grid.get());

  sim->setParameters(pList);

  sim->setup(); 

  if (sim->hasPEC()) {
    std::cout << "Saving PEC shape function...\n";
    shapes[pecShapeIndx]->save(sim->getGrid(), "pec");

#if 0
    MxShape<DIM> * sh;
    char name[200];
    for (int i = 0; i < shapes[pecShapeIndx]->getNumSubShapes(); i++) {
      sprintf(name, "pecSubShape%.2i", i);
      sh = shapes[pecShapeIndx]->getSubShape(i);
      sh->save(sim->getGrid(), name, name);
      delete sh;
    }
#endif
  }

  //solver = Teuchos::rcp(new MxSolver<DIM, MxComplex>(&sim, pList));
  //solver->setModeNames(eigModeNames);

  // solve the problem! This function also checks the eigensolutions
  // and saves the frequencies and fields to h5 files
  //solver->solve();

  bool isComplex = pList.get("is complex", false);
  if (isComplex) {
    MxSolver<DIM, MxComplex> solver(sim, pList);

    // solve the problem! This function also checks the eigensolutions
    // and saves the frequencies and fields to h5 files
    solver.solve();
  }
  else {
    MxSolver<DIM, double> solver(sim, pList);
    solver.solve();
  }



#if 0
  for (int i = 0; i < simXML.numChildren(); ++i)
    //parseNode(simXML.getChild(i));
    postParse(simXML.getChild(i));
#endif

  //if (pointClouds.size() != 0) {
  //  for (size_t i = 0; i < pointClouds.size(); ++i) {
  //    solver->saveFieldValsAtPoints(pointClouds[i]);
  //  }
  //}

  //Teuchos::RCP<Epetra_MultiVector> sol = solver->
  
  //setGeometry();

  //solve();

}

template<size_t DIM>
MxProblem<DIM>::~MxProblem() {
  for (size_t i = 0; i < shapes.size(); ++i)
    delete shapes[i];
}

template<size_t DIM>
void MxProblem<DIM>::preParse(const Teuchos::XMLObject & node) {
  std::string tag = node.getTag();

  if (shapeTags.find(tag) != shapeTags.end())
    getShape(node);
  else if (tag == "Eigensolver")
    getEigensolver(node);
  else if (tag == "Grid")
    getGrid(node);
  else if (tag == "Dielectric")
    getDielectric(node);
  else if (tag == "Mu")
    getMu(node);
  else if (tag == "PML")
    getPML(node);
  else if (tag == "PEC")
    getPEC(node);
  else if (tag == "BoundaryConditions")
    getBCs(node);
  else
    std::cout << "MxProblem<" << DIM << ">::parseNode: Did not understand XML tag: '" << tag << "'. Ignoring.\n";
}



template<size_t DIM>
void MxProblem<DIM>::parseNode(const Teuchos::XMLObject & node) {
  std::string tag = node.getTag();

  if (shapeTags.find(tag) != shapeTags.end())
    getShape(node);
  else if (tag == "Eigensolver")
    getEigensolver(node);
  else if (tag == "Grid")
    getGrid(node);
  else if (tag == "Dielectric")
    getDielectric(node);
  else if (tag == "Mu")
    getMu(node);
  else if (tag == "PEC")
    getPEC(node);
  else if (tag == "BoundaryConditions")
    getBCs(node);
  //else if (tag == "Output")
  //  getOutput(node);
  else
    std::cout << "MxProblem<" << DIM << ">::parseNode: Did not understand XML tag: '" << tag << "'. Ignoring.\n";

}

template<size_t DIM>
size_t MxProblem<DIM>::getShape(const Teuchos::XMLObject & node) {
  std::string tag = node.getTag();
  MxShape<DIM> * shape = 0;
  if (tag == "Ellipsoid")
    shape = new MxEllipsoid<DIM>(node);
  else if (tag == "Sphere")
    shape = new MxSphere<DIM>(node);
  else if (tag == "Cylinder")
    shape = new MxCylinder<DIM>(node);
  else if (tag == "Torus")
    shape = new MxTorus<DIM>(node);
  else if (tag == "Cone")
    shape = new MxCone<DIM>(node);
  else if (tag == "Slab")
    shape = new MxSlab<DIM>(node);
  else if (tag == "HalfSpace")
    shape = new MxHalfSpace<DIM>(node);
  else if (tag == "ShapeSubtract") {
    // node should have two children
    int n = node.numChildren();
    if (n != 2)
      std::cout << "'ShapeSubtract' should have exactly 2 shape subobjects. Found " << n << ".\n";
    else {
      MxShapeSubtract<DIM> * tmp = new MxShapeSubtract<DIM>();

      size_t baseIndx = getShape(node.getChild(0));
      size_t subtIndx = getShape(node.getChild(1));

      tmp->setBaseShape(shapes[baseIndx]);
      tmp->subtractShape(shapes[subtIndx]);

      shape = tmp;
      shape->transforms(node);
    }
  }
  else if (tag == "ShapeIntersection") {
    MxShapeIntersection<DIM> * tmp = new MxShapeIntersection<DIM>();
    size_t indx;
    for (int i = 0; i < node.numChildren(); ++i) {
      indx = getShape(node.getChild(i));
      tmp->add(shapes[indx]);
    }
    shape = tmp;
    shape->transforms(node);
  }
  else if (tag == "ShapeUnion") {
    MxShapeUnion<DIM> * tmp = new MxShapeUnion<DIM>();
    size_t indx;
    for (int i = 0; i < node.numChildren(); ++i) {
      indx = getShape(node.getChild(i));
      tmp->add(shapes[indx]);
    }
    shape = tmp;
    shape->transforms(node);
  }
  else if (tag == "ShapeMirror") {
    size_t indx;
    MxDimVector<double, DIM> n, p;
    n.strFill(MxUtil::XML::getAttr("normal", node));
    p.strFill(MxUtil::XML::getAttr("point in plane", node));
    if (node.numChildren() != 1) {
      std::cout << "ShapeMirror XML block must contain one shape block.\n";
      node.print(std::cout, 1);
      exit(0);
    }
    else {
      indx = getShape(node.getChild(0));
    }
    shape = new MxShapeMirror<DIM>(shapes[indx], n, p);
    // apply generic shape transformations (order matters in .mx file!)
    shape->transforms(node);
  }
  else if (tag == "ShapeRepeat") {
    size_t indx;
    MxDimVector<double, DIM> dir, o;
    double step;
    int numPos, numNeg;
    o.strFill(MxUtil::XML::getAttr("origin", node));
    dir.strFill(MxUtil::XML::getAttr("direction", node));
    step = atof(MxUtil::XML::getAttr("step", node).c_str());
    numPos = atoi(MxUtil::XML::getAttr("pos steps", node).c_str());
    numNeg = atoi(MxUtil::XML::getAttr("neg steps", node).c_str());
    if (node.numChildren() != 1) {
      std::cout << "ShapeRepeat XML block must contain one shape block.\n";
      node.print(std::cout, 1);
      exit(0);
    }
    else {
      indx = getShape(node.getChild(0));
    }
    shape = new MxShapeRepeat<DIM>(shapes[indx], o, dir, step, numPos, numNeg);
    // apply generic shape transformations (order matters in .mx file!)
    shape->transforms(node);
  }
  else
    std::cout << "MxProblem<" << DIM << ">::getShape: Did not understand shape tag, '" << tag << "'. Skipping object.\n";

  if (shape) {
    shapes.push_back(shape);
    return shapes.size() - 1;
  }
  else
    return -1;

}

template<size_t DIM>
void MxProblem<DIM>::getLinearSolver(const Teuchos::XMLObject & node) {
  if (pid == 0) std::cout << "Setting linear solver parameters...\n";

  for (int i = 0; i < numLinSolverParams; ++i)
    addParam(node, "linear solver : ", linSolverParams[i]);

  // search for a preconditioner
  Teuchos::XMLObject prec;
  int res = node.findFirstChild("Preconditioner");
  if (res != -1) {
    prec = node.getChild(res);

    // get type
    std::string type = MxUtil::XML::getAttr("type", prec);
    pList.set("linear solver : prec type", type);
    if (type == "amg") {
      if (pid == 0) std::cout << "  Found AMG preconditioner\n";
      for (int i = 0; i < numAMGParams; ++i)
        addParam(prec, "linear solver : amg : ", amgParams[i]);
    }
    else if (type == "ilu") {
      if (pid == 0) std::cout << "  Found ILU preconditioner\n";
      for (int i = 0; i < numILUParams; ++i)
        addParam(prec, "linear solver : ilu : ", iluParams[i]);
    }
    else {
      if (pid == 0) std::cout << "Did not understand Preconditioner type: '"
        << type << "'.\n\nOffending XML block:\n" << prec;
      exit(EXIT_FAILURE);
    }
  }
  

}

template<size_t DIM>
void MxProblem<DIM>::getEigensolver(const Teuchos::XMLObject & node) {
  if (pid == 0) std::cout << "Setting eigensolver parameters...\n";

  for (int i = 0; i < numEigSolverParams; ++i)
    addParam(node, "eigensolver : ", eigSolverParams[i]);

  bool invert = pList.get("eigensolver : invert", true);
  if (invert) {
    // search for LinearSolver block
    Teuchos::XMLObject linSolver;
    int res = node.findFirstChild("LinearSolver");
    if (res != -1)
      getLinearSolver(node.getChild(res));
    else {
      if (pid == 0)
      std::cout << "In eigensolver block, invert = true but no child <LinearSolver> block was given.\n";
      exit(EXIT_FAILURE);
    }
  }

  eigSet = true;
}

template<size_t DIM>
void MxProblem<DIM>::getGrid(const Teuchos::XMLObject & node) {
  MxDimVector<int, DIM> n;
  MxDimVector<double, DIM> l, o;

  n.strFill(MxUtil::XML::getAttr("resolution", node, "[16, 16, 16]"));
  l.strFill(MxUtil::XML::getAttr("dimensions", node, "[1.0, 1.0, 1.0]"));
  o.strFill(MxUtil::XML::getAttr("origin", node, "[0.0, 0.0, 0.0]"));

  grid = Teuchos::rcp(new MxGrid<DIM>(o, n, l, comm));
  grid->print();

  gridSet = true;
}

// a PEC object should contain only one Shape child object
template<size_t DIM>
void MxProblem<DIM>::getPEC(const Teuchos::XMLObject & node) {
  if (node.numChildren() > 1) {
    std::cout << "PEC object should contain no more than 1 Shape object";
    exit(0);
  }

  double dmFrac = atof(MxUtil::XML::getAttr("D-M frac", node, "0.0").c_str());
   
  pecShapeIndx = getShape(node.getChild(0));

  sim->setPEC(shapes[pecShapeIndx], dmFrac);

#if 0
  std::cout << "Number of sub-shapes in PEC: " << shapes[pecShapeIndx]->getNumSubShapes() << "\n";
  MxShape<DIM> * sh;
  for (int i = 0; i < shapes[pecShapeIndx]->getNumSubShapes(); ++i) {
    sh = shapes[pecShapeIndx]->getSubShape(i);
    std::cout << "  " << i << ": " << sh->getName() << "\n";
    delete sh;
  }
#endif
}


// should contain one shape block and some epsilon attributes
template<size_t DIM>
void MxProblem<DIM>::getDielectric(const Teuchos::XMLObject & node) {
  if (node.numChildren() > 1) {
    std::cout << "Dielectric object should contain no more than 1 Shape object";
    exit(0);
  }

  // get the name of this dielectric object (defaults to 'dielectric##')
  std::string name = MxUtil::XML::getAttr("name", node,
      std::string("dielectric") + MxUtil::Strings::typeToStr(sim->numDielectrics()));
  if (pid == 0) std::cout << "Parsed dielectric object '" << name << "'\n";
  MxDielectric<DIM> * diel = new MxDielectric<DIM>(name);

  size_t indx = getShape(node.getChild(0));
  dielShapeIndices.push_back(indx);

  MxDimVector<double, 3> diag, offDiag;
  //std::cout << "finding epses\n";
  diag.strFill(MxUtil::XML::getAttr("eps diag", node, "[1.0, 1.0, 1.0]"));
  offDiag.strFill(MxUtil::XML::getAttr("eps off diag", node, "[0.0, 0.0, 0.0]"));
  double losstan = atof(MxUtil::XML::getAttr("loss tangent", node, "0.0").c_str());

  if (losstan != 0)
    pList.set("is complex", true);

  if (pid == 0) {
    std::cout << "  Diagonal epsilon: " << diag;
    std::cout << "  Off-diagonal epsilon: " << offDiag;
    std::cout << "  Loss tangent: " << losstan << "\n";
    std::cout << "  Shape: " << shapes[indx]->getName() << "\n\n";
  }

  MxDimMatrix<MxComplex, 3> eps;
  for (size_t i = 0; i < 3; ++i)
    eps(i, i) = MxComplex(diag[i], diag[i] * losstan);
  eps(0, 1) = eps(1, 0) = MxComplex(offDiag[2], 0);
  eps(0, 2) = eps(2, 0) = MxComplex(offDiag[1], 0);
  eps(1, 2) = eps(2, 1) = MxComplex(offDiag[0], 0);

  diel->setEps(shapes[indx], eps);

  sim->addDielectric(diel);
}

// should contain one shape block and some epsilon attributes
template<size_t DIM>
void MxProblem<DIM>::getMu(const Teuchos::XMLObject & node) {
  if (node.numChildren() > 1) {
    std::cout << "Mu object should contain no more than 1 Shape object";
    exit(0);
  }

  // get the name of this dielectric object (defaults to 'dielectric##')
  std::string name = MxUtil::XML::getAttr("name", node,
      std::string("mu") + MxUtil::Strings::typeToStr(sim->numMus()));
  if (pid == 0) std::cout << "Parsed mu object '" << name << "'\n";
  MxMu<DIM> * muObj = new MxMu<DIM>(name);

  size_t indx = getShape(node.getChild(0));
  muShapeIndices.push_back(indx);

  MxDimVector<double, 3> diag, offDiag;
  //std::cout << "finding epses\n";
  diag.strFill(MxUtil::XML::getAttr("mu diag", node, "[1.0, 1.0, 1.0]"));
  offDiag.strFill(MxUtil::XML::getAttr("mu off diag", node, "[0.0, 0.0, 0.0]"));
  double losstan = atof(MxUtil::XML::getAttr("loss tangent", node, "0.0").c_str());

  if (losstan != 0)
    pList.set("is complex", true);

  if (pid == 0) {
    std::cout << "  Diagonal mu: " << diag;
    std::cout << "  Off-diagonal mu: " << offDiag;
    std::cout << "  Loss tangent: " << losstan << "\n";
    std::cout << "  Shape: " << shapes[indx]->getName() << "\n\n";
  }

  MxDimMatrix<MxComplex, 3> mu;
  for (size_t i = 0; i < 3; ++i)
    mu(i, i) = MxComplex(diag[i], diag[i] * losstan);
  mu(0, 1) = mu(1, 0) = MxComplex(offDiag[2], 0);
  mu(0, 2) = mu(2, 0) = MxComplex(offDiag[1], 0);
  mu(1, 2) = mu(2, 1) = MxComplex(offDiag[0], 0);

  muObj->setMu(shapes[indx], mu);

  sim->addMu(muObj);
}

template<size_t DIM>
void MxProblem<DIM>::getPML(const Teuchos::XMLObject & node) {
  if (node.numChildren() > 1) {
    if (pid == 0) std::cout << "PML object should contain no more than 1 Shape object";
    exit(0);
  }
  if (node.numChildren() == 0) {
    if (pid == 0) std::cout << "PML object needs to contain exactly 1 Shape object";
    exit(0);
  }

  // get the name of this PML object (defaults to 'PML##')
  std::string name = MxUtil::XML::getAttr("name", node,
      std::string("PML") + MxUtil::Strings::typeToStr(sim->numPMLs()));
  MxPML<DIM> * pml = new MxPML<DIM>(name);

  // get the shape of this PML object
  size_t indx = getShape(node.getChild(0));
  pmlShapeIndices.push_back(indx);
  pml->setShape(shapes[indx]);

  // get maximum alpha
  double alphaMax = atof(MxUtil::XML::getAttr("max alpha", node, "1.0").c_str());
  pml->setAlphaMax(alphaMax);

  // get profile exponent
  double exponent = atof(MxUtil::XML::getAttr("exponent", node, "2.0").c_str());
  pml->setExponent(exponent);

  // get profile endpoints (required attributes)
  MxDimVector<double, DIM> profBegin, profEnd;
  profBegin.strFill(MxUtil::XML::getAttr("profile begin", node));
  profEnd.strFill(MxUtil::XML::getAttr("profile end", node));
  pml->setProfileEndpoints(profBegin, profEnd);

  // finally add it the the simulation
  sim->addPML(pml);

  // if any pmls are added, the simulation needs complex fields
  pList.set("is complex", true);
}


template<size_t DIM>
void MxProblem<DIM>::getBCs(const Teuchos::XMLObject & node) {
  MxDimVector<std::string, DIM> lowerBCs(""), upperBCs("");
  MxDimVector<double, DIM> phaseShifts(0);
  lowerBCs.strFill(MxUtil::XML::getAttr("lower bcs", node, "[periodic, periodic, periodic]"), "");
  upperBCs.strFill(MxUtil::XML::getAttr("upper bcs", node, "[periodic, periodic, periodic]"), "");
  phaseShifts.strFill(MxUtil::XML::getAttr("phase shifts", node, "[0.0, 0.0, 0.0]"), 0);

  MxDimVector<MxBCType, DIM> typeLBCs(PERIODIC), typeUBCs(PERIODIC);
  for (size_t i = 0; i < DIM; ++i) {
    // lower bcs
    if (lowerBCs[i] == "periodic")
      typeLBCs[i] = PERIODIC;
    else if (lowerBCs[i] == "pec")
      typeLBCs[i] = PEC;
    else if (lowerBCs[i] == "pmc")
      typeLBCs[i] = PMC;
    else {
      std::cout << "Invalid lower boundary condition: '" << lowerBCs[i] << "'.\n";
      exit(0);
    }

    // upper bcs
    if (upperBCs[i] == "periodic")
      typeUBCs[i] = PERIODIC;
    else if (upperBCs[i] == "pec")
      typeUBCs[i] = PEC;
    else if (upperBCs[i] == "pmc")
      typeUBCs[i] = PMC;
    else {
      std::cout << "Invalid upper boundary condition: '" << upperBCs[i] << "'.\n";
      exit(0);
    }
  }

  pList.set("upper bcs", typeUBCs);
  pList.set("lower bcs", typeLBCs);
  pList.set("phase shifts", phaseShifts);

  for (size_t i = 0; i < DIM; ++i)
    if (phaseShifts[i] != 0)
      pList.set("is complex", true);
}

#if 0
template<size_t DIM>
void MxProblem<DIM>::postParse(const Teuchos::XMLObject & node) {
  std::string tag = node.getTag();

  if (tag == "Output")
    getOutput(node);
  //else
  //  std::cout << "MxProblem<" << DIM << ">::parseNode: Did not understand XML tag: '" << tag << "'. Ignoring.\n";
}

template<size_t DIM>
void MxProblem<DIM>::getOutput(const Teuchos::XMLObject & node) {
  int numCh = node.numChildren();

  std::string tag;

  Teuchos::XMLObject child;
  for (int i = 0; i < numCh; ++i) {
    child = node.getChild(i);
    tag = child.getTag();
    if (tag == "FieldValues") {
      MxIO<DIM> io(&grid->getComm());

      Teuchos::RCP<Epetra_MultiVector> eevec(solver->getElecEigVecs());
      Teuchos::RCP<Epetra_MultiVector> bevec(solver->getMagEigVecs());

      std::string name = MxUtil::XML::getAttr("name", child);
      name += std::string("_eigenvectors");
      std::string kind = MxUtil::XML::getAttr("kind", child);
      if (kind == "points") {
        MxPointCloud<DIM> ptCld(child.getChild(0), grid.get());

        io.save(ptCld, *eevec, sim->getField("efield"), name);
        io.save(ptCld, *bevec, sim->getField("bfield"), name);
      }
      else if (kind == "grid") {
        MxGrid<DIM> tmpGrid(child.getChild(0), &grid->getComm());

        io.save(tmpGrid, *eevec, sim->getField("efield"), name);
        io.save(tmpGrid, *bevec, sim->getField("bfield"), name);
      }
      else {
        std::cout << "In FieldValues output, kind, '" << kind << ",' not understood.\n";
      }
    }
    else {
      std::cout << "While parsing <Output> block, child block with tag, '" << tag << ",' not understood.\n";
    }
  }
}

#endif

template<size_t DIM>
void MxProblem<DIM>::solve() {}


template class MxProblem<1>;
template class MxProblem<2>;
template class MxProblem<3>;
