
// This is for Epetra_ConfigDefs.h
//#define HAVE_CONFIG_H
//
//#include "MxConfig.h"
#include "Epetra_ConfigDefs.h"

// include mpi must be first for nersc machines
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <cmath>

#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxUtil.hpp"
#include "MxGrid.h"
#include "MxEllipsoid.hpp"
#include "MxShapeSubtract.hpp"
#include "MxSegment.h"
//#include "MxParallelogram.h"
#include "MxCartBox.hpp"
#include "MxCartRect.hpp"
//#include "MxCrabCav.h"

#include "MxDielectric.hpp"

#include "MxEMSim.h"
#include "MxSolver.h"
#include "MxProblem.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_XMLObject.hpp"

#include "Epetra_CrsMatrix.h"

int main(int argc, char* argv[]) {

#ifdef HAVE_MPI
  MPI::Init(argc, argv);
#endif
  RCP<MxComm> myComm = rcp(new MxComm());

#if 0
#ifdef HAVE_MPI
  MPI::Init(argc, argv);
  //MPI_Init(argc, argv);
  Epetra_MpiComm myComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm myComm;
#endif
#endif

// input file method
#if 1

  std::string inFile;

  Teuchos::CommandLineProcessor cmdp(false, true);
  cmdp.setOption("infile", &inFile, "XML format input file.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  if (inFile == "") {
    std::cout << "Please specify an input file using --infile=your_file.mx\n";
    exit(0);
  }

  // now read the input file with trilinos XML reader
  Teuchos::XMLObject xmlObj(Teuchos::FileInputSource(inFile).getObject());

  // get simulation dimension
  int dim = 0;
  // std::cout << "xmlObject = " << xmlObj << std::endl;
  try {
    dim = atoi(MxUtil::XML::getAttr("dim", xmlObj).c_str());
  } catch(std::exception e) {
    std::cout << "xmlObject = " << xmlObj << std::endl;
  }
  if (dim < 1 or dim > 3) {
    std::cout << "Simulation dimension invalid or not given, using 3D.\n";
    dim = 3;
  }

  // get simulation type
  std::string domain = MxUtil::XML::getAttr("domain", xmlObj).c_str();
  if (domain != "frequency" and domain != "time") {
    std::cout << "Simulation domain invalid or not given, using frequency-domain.\n";
    domain = "frequency";
  }

  // create problem

  MxProblem<1> * prob1d;
  MxProblem<2> * prob2d;
  MxProblem<3> * prob3d;
  switch (dim) {
    case 1:
      prob1d = new MxProblem<1>(xmlObj, myComm);
      prob1d->solve();
      delete prob1d;
      break;
    case 2:
      prob2d = new MxProblem<2>(xmlObj, myComm);
      prob2d->solve();
      delete prob2d;
      break;
    case 3:
      prob3d = new MxProblem<3>(xmlObj, myComm);
      prob3d->solve();
      delete prob3d;
      break;
  }
#endif


#if 0

  // epetra stuff test
  MxMap map(10, 0, myComm);
  Epetra_CrsMatrix mat(Copy, map, 0);
  int ind = 2;
  double val = 0;
  mat.InsertGlobalValues(1, 1, &val, &ind);
  ind = 3;
  val = 4;
  mat.InsertGlobalValues(1, 1, &val, &ind);
  mat.FillComplete(map, map);

  Epetra_Vector myvec(map);
  myvec.Random();

  std::cout << myvec;
  mat.Apply(myvec, myvec);
  std::cout << myvec;

  Epetra_CrsMatrix copy(mat);

  std::cout << mat;
  MxUtil::Epetra::stripZeros(mat);
  std::cout << mat;

  //throw 1;
#endif

  typedef MxDimVector<double, 3> vecd3;
  typedef MxDimVector<int, 3> veci3;

  vecd3 midPt(0);

#if 0
  //std::cout << "Crab cavity setup:\n";

  int crabNumCells = 4;
  double crabCellLen = 2.0 * 0.0192; //meters
  double crabCavRad = 0.04719;
  double crabIrisRad = 0.015;
  double crabCavRho = 0.0136;
  double crabIrisRho = 0.00331;

  int crabCellRes = 40;
  int padCells = 2;
  int cnx, cny, cnz;
  double clx, cly, clz;
  double cox, coy, coz;

  double crabDelta = crabCellLen / double(crabCellRes);

  cnz = crabNumCells * crabCellRes + 2 * padCells;
  clz = double(cnz) * crabDelta;
  coz = -0.5 * clz;

  cny = cnx = 2 * (int(ceil(crabCavRad / crabDelta)) + padCells);
  cly = clx = double(cnx) * crabDelta;
  coy = cox = -0.5 * clx;

  veci3 crabN; crabN[0] = cnx; crabN[1] = cny; crabN[2] = cnz;
  vecd3 crabL; crabL[0] = clx; crabL[1] = cly; crabL[2] = clz;
  vecd3 crabO; crabO[0] = cox; crabO[1] = coy; crabO[2] = coz;
  //crabN.print();
  //crabL.print();
  //crabO.print();

  MxGrid<3> crabGrid(crabO, crabN, crabL, &myComm);
  crabGrid.print();

  MxCrabCav crabCav(midPt, crabNumCells, crabCellLen, crabIrisRad, crabCavRad, crabIrisRho, crabCavRho);

  crabCav.save(crabGrid);

  Teuchos::ParameterList crabList;
  crabList.set("geo-mg : levels", 1);
  crabList.set("geo-mg : smoothers : sweeps", 5);
  crabList.set("amg : smoothers : sweeps", 1);
  crabList.set("amg : smoothers : type", "Chebyshev");
  crabList.set("eigensolver : output", 2);
  crabList.set("eigensolver : nev", 15);
  crabList.set("eigensolver : tol", 1.e-8);
  crabList.set("eigensolver : block size", 2);
  crabList.set("eigensolver : num blocks", 30);
  crabList.set("eigensolver : spectrum", "LM");
  crabList.set("wave operator : invert", true);
  crabList.set("wave operator : invert : tol", 1.e-10);
  //crabList.set("wave operator : invert : shift", 1000.0);
  crabList.set("wave operator : invert : max basis size", 40);

  MxEMSim<dim> crabSim;
  crabSim.setGrid(&crabGrid);
  crabSim.setPEC(&crabCav);
  //crabSim.setGrid(&sphGrid);
  //crabSim.setPEC(&ell);
  crabSim.setParameters(crabList);
  crabSim.setup();

  MxSolver<dim> * solver;
  solver = new MxSolver<dim>(&crabSim, crabList);
  solver->solve();

  delete solver;

  //return 1;
#endif

// optimized phc cavity
#if 0
  double rodRad = 0.003175; // meters
  const int numRods = 24;
  double rodx[numRods] = {0.0158406582694, 0.0551748491968, 0.0209567636489,
                          0.0384658321918, 0.00792032913471, 0.0338604938991,
                          0.00477355412058, 0.00485955186622, -0.00792032913471,
                          -0.0213143552977, -0.0161832095283, -0.0336062803256,
                          -0.0158406582694, -0.0551748491968, -0.0209567636489,
                          -0.0384658321918, -0.00792032913471, -0.0338604938991,
                          -0.00477355412058, -0.00485955186622, 0.00792032913471,
                          0.0213143552977, 0.0161832095283, 0.0336062803256};
  double rody[numRods] = {0.0, -0.00724351649877, 0.006587367621,
                    0.0165969314144, 0.013718412474, 0.044161062805,
                    0.0214427735115, 0.041610853563, 0.013718412474,
                    0.0514045793038, 0.0148554058905, 0.0250139221487,
                    1.9399211446e-18, 0.00724351649877, -0.006587367621,
                    -0.0165969314144, -0.013718412474, -0.044161062805,
                    -0.0214427735115, -0.041610853563, -0.013718412474,
                    -0.0514045793038, -0.0148554058905, -0.0250139221487};

  std::vector<MxShape<3> *> rods;
  MxShapeUnion<3> rodsShape;
  vecd3 rodPos;
  vecd3 zhat(0); zhat[2] = 1.0;
  for (int i = 0; i < numRods; i++) {
    rodPos[0] = rodx[i];
    rodPos[1] = rody[i];
    rodPos[2] = 0.0;
    rods.push_back(new MxCylinder(rodPos, zhat, rodRad));
    rodsShape.add(rods[i]);
  }

  MxDimMatrix<double, 3> sapphEps(0);
  sapphEps(0, 0) = 9.3;
  sapphEps(1, 1) = 9.3;
  sapphEps(2, 2) = 11.5;

  MxDielectric<3> phcDiel;
  phcDiel.add(&rodsShape, sapphEps);

  // conducting cavity
  double cavLen = 0.019624116824498831;
  double cavRad = 0.1;
  MxCylinder cavCyl(0, zhat, cavRad);
  MxSlab<3> cavCaps(0, zhat, cavLen);
  MxShapeIntersection<3> phcCav;
  phcCav.add(&cavCyl);
  phcCav.add(&cavCaps);

  // setup grid
  int rodDiaCells = 6;
  int pad = 2;
  double delta = 2.0 * rodRad / double(rodDiaCells);

  veci3 phcN;
  phcN[0] = phcN[1] = int(2.0 * cavRad / delta) + 2 * pad;
  phcN[2] = int(cavLen / delta) + 2 * pad;

  vecd3 phcL;
  phcL[0] = phcL[1] = delta * double(phcN[0]);
  phcL[2] = delta * double(phcN[2]);

  vecd3 phcO;
  phcO[0] = phcO[1] = -0.5 * phcL[0];
  phcO[2] = -0.5 * phcL[2];

  MxGrid<3> phcGrid(phcO, phcN, phcL, &myComm);
  phcGrid.print();

  Teuchos::ParameterList phcList;
  phcList.set("geo-mg : levels", 1);
  phcList.set("geo-mg : smoothers : sweeps", 5);
  phcList.set("eigensolver : output", 2);
  phcList.set("eigensolver : nev", 15);
  phcList.set("eigensolver : tol", 1.e-8);
  phcList.set("eigensolver : block size", 1);
  phcList.set("eigensolver : num blocks", 30);
  phcList.set("eigensolver : spectrum", "LM");
  phcList.set("wave operator : invert", true);
  phcList.set("wave operator : invert : tol", 1.e-8);
  //phcList.set("wave operator : invert : shift", 1000.0);
  phcList.set("wave operator : invert : max basis size", 40);

  MxEMSim<dim> phcSim;
  phcSim.setGrid(&phcGrid);
  //phcSim.setPEC(&phcCav);
  phcSim.setDielectric(&phcDiel);
  phcSim.setParameters(phcList);
  phcSim.setup();

  MxSolver<dim> * solver;
  solver = new MxSolver<dim>(&phcSim, phcList);
  solver->solve();

  delete solver;

  for (int i = 0; i < numRods; i++)
    delete rods[i];

#endif

#if 0
  double sphR = 0.37;
  int sphN = 64;
  MxEllipsoid ell(0.0, sphR);

  MxGrid<3> sphGrid(-0.5, sphN, 1.0, &myComm);
  sphGrid.print();

  MxDimMatrix<double, 3> rotSapphEps(0);
  rotSapphEps(0, 0) = 10.225;
  rotSapphEps(1, 1) = 10.225;
  rotSapphEps(2, 2) = 9.95;
  rotSapphEps(0, 1) = rotSapphEps(1, 0) = -0.825;
  rotSapphEps(0, 2) = rotSapphEps(2, 0) = -0.67360967926537398;
  rotSapphEps(1, 2) = rotSapphEps(2, 1) = 0.67360967926537398;

  MxDielectric<3> phcDiel;
  phcDiel.add(&ell, rotSapphEps);

  vecd3 ell2Loc(0); ell2Loc[0] = 0.6;
  vecd3 ell3Loc(0); ell3Loc[0] = 0.3; ell3Loc[2] = 0.3;
  MxEllipsoid ell2(ell2Loc, sphR);
  MxEllipsoid ell3(ell3Loc, sphR);

  MxShapeUnion<3> shUnion;
  shUnion.add(&ell);
  shUnion.add(&ell2);
  shUnion.add(&ell3);
  //shUnion.save(sphGrid);

  MxShapeIntersection<3> shInt;
  shInt.add(&ell);
  shInt.add(&ell2);
  shInt.add(&ell3);
  //shInt.save(sphGrid);

  MxShapeSubtract<3> shSub;
  shSub.setBaseShape(&ell);
  shSub.subtractShape(&ell2);
  shSub.subtractShape(&ell3);
  //shSub.save(sphGrid);

  MxDielectric<3> dielEll;
  MxDimMatrix<double, 3> epsEll(vecd3(10.0)); // isotropic eps = 10
  dielEll.add(&ell, epsEll);

  Teuchos::ParameterList sphList;
  sphList.set("geo-mg : levels", 1);
  sphList.set("geo-mg : smoothers : sweeps", 4);
  sphList.set("eigensolver : output", 2);
  sphList.set("eigensolver : nev", 12);
  sphList.set("eigensolver : tol", 1.e-8);
  sphList.set("eigensolver : block size", 1);
  sphList.set("eigensolver : num blocks", 30);
  sphList.set("eigensolver : spectrum", "LM");
  sphList.set("wave operator : invert", true);
  sphList.set("wave operator : invert : tol", 1.e-8);
  //sphList.set("wave operator : invert : shift", -0.1);
  sphList.set("wave operator : invert : shift", 1.0);
  sphList.set("wave operator : invert : max basis size", 40);

  MxEMSim<dim> sphSim;
  sphSim.setGrid(&sphGrid);
  //sphSim.setDielectric(&dielEll);
  sphSim.setDielectric(&phcDiel);
  //sphSim.setPEC(&sphCav);
  //sphSim.setPEC(&ell);
  sphSim.setParameters(sphList);
  sphSim.setup();

  MxSolver<dim> * solver;
  solver = new MxSolver<dim>(&sphSim, sphList);
  solver->solve();

  delete solver;

#endif


#ifdef HAVE_MPI
  MPI::Finalize();
  //MPI_Finalize();
#endif

  return 0;

}
