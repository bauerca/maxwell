
//#include "Ifpack_ConfigDefs.h"

#include "MxGeoMultigridPrec.h"

#include <ctime>

#include "MxUtil.hpp"

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Ifpack.h"
#include "Ifpack_Amesos.h"
#include "Ifpack_Chebyshev.h"

MxGeoMultigridPrec::MxGeoMultigridPrec(Teuchos::ParameterList aPList) :
pList(aPList), totalTime(0.0), numApplyInverse(0) {}

MxGeoMultigridPrec::~MxGeoMultigridPrec() {
  std::cout << "============= MxGeoMultigridPrec =============\n"
            << "  ApplyInverse called " << numApplyInverse << "times\n"
            << "  Total time taken: " << totalTime << " sec\n"
            << "  Time per ApplyInverse call: " << totalTime / double(numApplyInverse) << " sec\n"
            << "==============================================\n";
}

void MxGeoMultigridPrec::setOperators(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs) {
  ops = matPtrs;
  opsSet = true;
  levels = ops.size();
}

void MxGeoMultigridPrec::setCoarseners(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs) {
  setLHSCoarseners(matPtrs);
  setRHSCoarseners(matPtrs);
}

void MxGeoMultigridPrec::setRHSCoarseners(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs) {
  if (!opsSet) {
    std::cout << "MxGeoMultigridPrec::setRHSCoarseners: must set operators first.\n";
    throw 1;
  }
  if (matPtrs.size() != ops.size() - 1) {
    std::cout << "MxGeoMultigridPrec::setRHSCoarseners: number of coarseners: " << matPtrs.size()
              << ", number of operators: " << levels << ".\n";
    throw 1;
  }
  rhsCoarseners = matPtrs;
}

void MxGeoMultigridPrec::setLHSCoarseners(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs) {
  if (!opsSet) {
    std::cout << "MxGeoMultigridPrec::setLHSCoarseners: must set operators first.\n";
    throw 1;
  }
  if (matPtrs.size() != ops.size() - 1) {
    std::cout << "MxGeoMultigridPrec::setLHSCoarseners: number of coarseners: " << matPtrs.size()
              << ", number of operators: " << levels << ".\n";
    throw 1;
  }
  lhsCoarseners = matPtrs;
}

void MxGeoMultigridPrec::setRefiners(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs) {
  setLHSRefiners(matPtrs);
  setRHSRefiners(matPtrs);
}

void MxGeoMultigridPrec::setLHSRefiners(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs) {
  if (!opsSet) {
    std::cout << "MxGeoMultigridPrec::setLHSRefiners: must set operators first.\n";
    throw 1;
  }
  if (matPtrs.size() != ops.size() - 1) {
    std::cout << "MxGeoMultigridPrec::setLHSRefiners: number of refiners: " << matPtrs.size()
              << ", number of operators: " << levels << ".\n";
    throw 1;
  }
  lhsRefiners = matPtrs;
}

void MxGeoMultigridPrec::setRHSRefiners(std::vector<Teuchos::RCP<Epetra_CrsMatrix> > matPtrs) {
  if (!opsSet) {
    std::cout << "MxGeoMultigridPrec::setRHSRefiners: must set operators first.\n";
    throw 1;
  }
  if (matPtrs.size() != ops.size() - 1) {
    std::cout << "MxGeoMultigridPrec::setRHSRefiners: number of refiners: " << matPtrs.size()
              << ", number of operators: " << levels << ".\n";
    throw 1;
  }
  rhsRefiners = matPtrs;
}


void MxGeoMultigridPrec::setup() {

  std::string smoothType = pList.get("linear solver : smoother type", "Gauss-Seidel");
  std::string coarseSmoothType = pList.get("linear solver : coarse smoother", "Amesos");
  int smoothSweeps = pList.get("linear solver : smoother sweeps", 1);
  int overlap = 1;
  
  levels = pList.get("linear solver : levels", 1);
  cycles = pList.get("linear solver : cycles", 1);
  output = pList.get("linear solver : output", 0);
  rcf = pList.get("linear solver : remove const field", false);

  // Gauss-Seidel preconditioner options
  Teuchos::ParameterList GSList;
  GSList.set("relaxation: type", "Gauss-Seidel");
  //GSList.set("relaxation: type", "Jacobi");
  GSList.set("relaxation: sweeps", smoothSweeps);
  GSList.set("relaxation: damping factor", 1.);
  GSList.set("relaxation: zero starting solution", false);

  // Chebyshev preconditioner options
  Teuchos::ParameterList ChebList;
  ChebList.set("chebyshev: degree", smoothSweeps);
  ChebList.set("chebyshev: zero starting solution", false);
  ChebList.set("chebyshev: ratio eigenvalue", 30.);

  // Amesos direct solver parameters
  Teuchos::ParameterList AmesosList;
  AmesosList.set("amesos: solver type", "Amesos_Klu");
  AmesosList.set("AddToDiag", 1.e-12);



  // smoother parameters for intermediate levels
  Teuchos::ParameterList vSmootherList;
  std::string vPrecType;
  if (smoothType == "Gauss-Seidel") {
    vPrecType = "point relaxation";
    vSmootherList = GSList;
  }
  else if (smoothType == "Chebyshev") {
    vPrecType = smoothType;
    vSmootherList = ChebList;
  }
  else {
    std::cout << "MxGeoMultigridPrec::setup(): unknown smoother type, '" << smoothType << "' for intermediate smoothers. Using gauss-seidel.\n";
    vSmootherList = GSList;
  }

  // smoother parameters for coarsest level
  Teuchos::ParameterList cSmootherList;
  std::string cPrecType;
  if (coarseSmoothType == "Gauss-Seidel") {
    cSmootherList = GSList;
    cPrecType = "point relaxation";
  }
  else if (coarseSmoothType == "Chebyshev") {
    cSmootherList = ChebList;
    cPrecType = coarseSmoothType;
  }
  else if (coarseSmoothType == "Amesos") {
    cSmootherList = AmesosList;
    cPrecType = coarseSmoothType;
  }
  else {
    std::cout << "MxGeoMultigridPrec::setup(): unknown smoother type, '" << coarseSmoothType << "', for coarsest smoother. Using Amesos_Klu.\n";
    cSmootherList = AmesosList;
    cPrecType = "Amesos";
  }

  Teuchos::ParameterList smootherList;
  std::string precType;

  Ifpack factory;
  
  double maxEig, minEig;

  // special case for when Amesos wants to alter our original operator
  if (levels == 1 and coarseSmoothType == "Amesos") {
    fineOpCopy = Teuchos::rcp(new Epetra_CrsMatrix(*ops[0])); //apparently Amesos changes the input matrix?
    smoothers.push_back(Teuchos::rcp(new Ifpack_Amesos(&*fineOpCopy)));
    smoothers[0]->SetParameters(AmesosList);
    smoothers[0]->Initialize();
    smoothers[0]->Compute();
  }
  else {
  //level 0 is original operator
    for (int level = 0; level < levels; ++level) {
      std::cout << "  Initializing multigrid level " << level << std::endl;

      if (level == levels - 1) {
        smootherList = cSmootherList;
        precType = cPrecType;
      }
      else {
        smootherList = vSmootherList;
        precType = vPrecType;
      }

      // always find max eigenvalue for at least the prolongator smoother,
      // and possibly the cheyshev smoother
      {
        Epetra_Vector diag(ops[level]->RangeMap());
        ops[level]->ExtractDiagonalCopy(diag);
        diag.Reciprocal(diag);

        std::cout << "    Inf norm: " << ops[level]->NormInf() << "\n";

        //diag.PutScalar(1.0);
        Ifpack_Chebyshev::CG(*ops[level], diag, 100, minEig, maxEig);
        std::cout << "unscaled operator: lambda max = " << maxEig 
                                   << ", lambda min: " << minEig << "\n";
        maxEigs.push_back(maxEig);

        if (precType == "Chebyshev") {
          smootherList.set("chebyshev: max eigenvalue", maxEig);
          smootherList.set("chebyshev: min eigenvalue", minEig);
        }
      }

      smoothers.push_back(Teuchos::rcp(factory.Create(precType, &*ops[level], overlap)));
      smoothers[level]->SetParameters(smootherList);
      smoothers[level]->Initialize();
      smoothers[level]->Compute();

      std::cout << *smoothers[level];

#if 0
      // check klu solver
      if (level == levels - 1) {
        std::cout << "Checking KLU solver with Ax = 0\n";
        Epetra_MultiVector zero(ops[level]->RowMap(), 1), work(ops[level]->RowMap(), 1);
        zero.PutScalar(0.0);
        work.Random();
        smoothers[level]->ApplyInverse(zero, work);
        std::cout << work;
      }
#endif


    }
  }
} 


int MxGeoMultigridPrec::vCycle(std::vector<Epetra_MultiVector*> & bvecs, 
                         std::vector<Epetra_MultiVector*> & xvecs,
                         int startlevel, int cycle) const {

  // output vars
  //double totalRes, bndryRes, bulkRes;
  //char name[200];

  // ptrs
  Epetra_MultiVector *b, *x;

  //workVecs2[0]->Random(); //initial guess
  // downstroke
  for (int i = startlevel; i < levels; ++i) {
    b = bvecs[i];
    x = xvecs[i];

    if (i != startlevel) 
      x->PutScalar(0.0);

#if 0
    if (output_) {
      spaces(i); std::cout << "downstroke pre smooth:\n";
      GetBoundaryBulkResiduals(*(ops_[i]), *x, *b, bndryRes, bulkRes, totalRes, i);
      //std::cout << "    residual before smooth: " << PhcGeoMGPrec::GetResidual(*(ops_[i]), *wv2, *wv1) << "\n";
      spaces(i); std::cout << "  total residual: " << totalRes << "\n";
      spaces(i); std::cout << "  bndry residual: " << bndryRes << "\n";
      spaces(i); std::cout << "  bulk  residual: " << bulkRes << "\n";

      // pre relaxation (or full solve if lowest level)
      //
      if (i==startlevel)
        sprintf(name, "geoMG-V%i-cycle%i-level%i-x-pre-down.h5", startlevel, cycle, i);
      else
        sprintf(name, "geoMG-V%i-cycle%i-level%i-e-pre-down.h5", startlevel, cycle, i);
      spaces(i); std::cout << "Saving " << name << "...";
      bfields_[i]->SaveToH5(name, comm_, x);
      if (i==startlevel)
        sprintf(name, "geoMG-V%i-cycle%i-level%i-b.h5", startlevel, cycle, i);
      else
        sprintf(name, "geoMG-V%i-cycle%i-level%i-r.h5", startlevel, cycle, i);
      spaces(i); std::cout << "Saving " << name << "...";
      bfields_[i]->SaveToH5(name, comm_, b);
    }
#endif

    //spaces(i); std::cout << "residual before smooth: " << getResidual(*ops[i], *x, *b) << "\n";
    smoothers[i]->ApplyInverse(*b, *x);
    //spaces(i); std::cout << "residual after smooth: " << getResidual(*ops[i], *x, *b) << "\n";

#if 0
    if (output_) {
      if (i==startlevel)
        sprintf(name, "geoMG-V%i-cycle%i-level%i-x-post-down.h5", startlevel, cycle, i);
      else
        sprintf(name, "geoMG-V%i-cycle%i-level%i-e-post-down.h5", startlevel, cycle, i);
      bfields_[i]->SaveToH5(name, comm_, x);

      spaces(i); std::cout << "downstroke post smooth:\n";
      GetBoundaryBulkResiduals(*(ops_[i]), *x, *b, bndryRes, bulkRes, totalRes, i);
      spaces(i); std::cout << "  total residual: " << totalRes << "\n";
      spaces(i); std::cout << "  bndry residual: " << bndryRes << "\n";
      spaces(i); std::cout << "  bulk  residual: " << bulkRes << "\n";
      //std::cout << "    residual after smooth: " << PhcGeoMGPrec::GetResidual(*(ops_[i]), *wv2, *wv1) << "\n";
    }
#endif

    // compute and coarsen residual (r = b - A\tilde{x})
    if (i < levels - 1) {
      Epetra_MultiVector r(*b);
      ops[i]->Apply(*x, r);
      r.Update(1., *b, -1.); //r=b-Ax

      //interpolators_[i]->SetUseTranspose(true);             //full-weighting coarsening
      //interpolators_[i]->Apply(worktmp, *(workVecs1[i+1])); //full-weighting coarsening
      //workVecs1[i+1]->Scale(0.125);                         //full-weighting coarsening

      // to coarsen: scale by DMAreas^{-1} to get true fields, 
      //             coarsen
      //             scale by DMAreas
      //worktmp.ReciprocalMultiply(1., *dmAreas_[i], worktmp, 0.);
      //coarseners[i]->Apply(r, *bvecs[i + 1]); //coarsened residual becomes new b for next level
      coarsen(i + 1, r, *bvecs[i + 1]);
      //SmoothInterpolation(i+1, *bvecs[i+1]); // b -> (I - op/rhomax)*b

      //workVecs1[i+1]->Multiply(1., *dmAreas_[i+1], *workVecs1[i+1], 0.);
    }
  }

  // upstroke
  for (int i = levels - 2; i >= startlevel; --i) {
    b = bvecs[i]; //holds b
    x = xvecs[i]; //holds x
    Epetra_MultiVector e(*b);
    // interpolate error field
    //interpolators_[i]->SetUseTranspose(false);
    //
    // to finen: scale by DMAreas to get fluxes, 
    //             interpolate fluxes
    //             scale by DMAreas^{-1}
    //workVecs2[i+1]->Multiply(1., *dmAreas_[i+1], *workVecs2[i+1], 0.);
    //refiners[i]->Apply(*xvecs[i+1], e); //vtmp holds interpolated e
    refine(i, *xvecs[i + 1], e); //vtmp holds interpolated e
    //vtmp.ReciprocalMultiply(1., *dmAreas_[i], vtmp, 0.);
    //SmoothInterpolation(i, e);  // e = (I - op/rhomax)*e

    // correct with interpolated error field
    x->Update(1., e, 1.); // x = x + e

    // post relaxation
    //std::cout << "    residual before smooth: " << PhcGeoMGPrec::GetResidual(*(ops_[i]), *wv2, *wv1) << "\n";
    // Print residuals before applying smoother
#if 0
    if (output_) {
      spaces(i); std::cout << "upstroke pre smooth:\n";
      GetBoundaryBulkResiduals(*(ops_[i]), *x, *b, bndryRes, bulkRes, totalRes, i);
      //std::cout << "    residual before smooth: " << PhcGeoMGPrec::GetResidual(*(ops_[i]), *wv2, *wv1) << "\n";
      spaces(i); std::cout << "  total residual: " << totalRes << "\n";
      spaces(i); std::cout << "  bndry residual: " << bndryRes << "\n";
      spaces(i); std::cout << "  bulk  residual: " << bulkRes << "\n";

      // save fields before applying smoother
      if (i==startlevel)
        sprintf(name, "geoMG-V%i-cycle%i-level%i-x-pre-up.h5", startlevel, cycle, i);
      else
        sprintf(name, "geoMG-V%i-cycle%i-level%i-e-pre-up.h5", startlevel, cycle, i);
      bfields_[i]->SaveToH5(name, comm_, x);
    }
#endif

    // apply smoother
    //spaces(i); std::cout << "residual before smooth: " << getResidual(*ops[i], *x, *b) << "\n";
    smoothers[i]->ApplyInverse(*b, *x);
    //spaces(i); std::cout << "residual after smooth: " << getResidual(*ops[i], *x, *b) << "\n";

#if 0
    if (output_) {
      // save fields after applying smoother
      if (i==startlevel)
        sprintf(name, "geoMG-V%i-cycle%i-level%i-x-post-up.h5", startlevel, cycle, i);
      else
        sprintf(name, "geoMG-V%i-cycle%i-level%i-e-post-up.h5", startlevel, cycle, i);
      bfields_[i]->SaveToH5(name, comm_, x);

      // print residuals after applying smoother
      spaces(i); std::cout << "upstroke post smooth:\n";
      GetBoundaryBulkResiduals(*(ops_[i]), *x, *b, bndryRes, bulkRes, totalRes, i);
      spaces(i); std::cout << "  total residual: " << totalRes << "\n";
      spaces(i); std::cout << "  bndry residual: " << bndryRes << "\n";
      spaces(i); std::cout << "  bulk  residual: " << bulkRes << "\n";
    }
#endif
  }

  return 0;
}

void MxGeoMultigridPrec::removeConstField(Epetra_MultiVector & x) const {
  //x.Comm().Barrier();
  Epetra_MultiVector ones(x);
  double dotProd, norm;
  ones.PutScalar(1);
  ones.Norm2(&norm);
  x.Dot(ones, &dotProd);
  //std::cout << "norm = " << norm << ", dotProd = " << dotProd << "\n";
  x.Update(-dotProd / norm / norm, ones, 1.);
  x.Dot(ones, &dotProd);
  //std::cout << "const vec part: " << dotProd << "\n";
}

double MxGeoMultigridPrec::constPart(const Epetra_MultiVector & x) const {
  //x.Comm().Barrier();
  Epetra_MultiVector ones(x);
  double dotProd;
  ones.PutScalar(1);
  x.Dot(ones, &dotProd);
  return dotProd;
}

// return  x = (I - op/rhomax)*x in x
void MxGeoMultigridPrec::smoothInterpolation(int level, Epetra_MultiVector & x) const {
  Epetra_MultiVector xtmp(x);

  // get matrix diagonal
  Epetra_Vector diag(ops[level]->RangeMap());
  ops[level]->ExtractDiagonalCopy(diag);
  //diag.Reciprocal(diag);

  ops[level]->Apply(x, xtmp);    // xtmp = op*x
  xtmp.ReciprocalMultiply(1.0, diag, xtmp, 0.0);
  //xtmp.Scale(1.3333 / maxEigs[level]); // xtmp = (op/rhomax)*x
  xtmp.Scale(0.5 * 1.3333); // xtmp = (op/rhomax)*x (rhomax is usually 2)
  x.Update(-1.0, xtmp, 1.0);        // x    = x - xtmp
}

void MxGeoMultigridPrec::refine(int targetLevel, const Epetra_MultiVector & coarse, Epetra_MultiVector & fine) const {
  lhsRefiners[targetLevel]->Apply(coarse, fine);
  if (rcf) MxUtil::Epetra::removeConstField(fine);

  //smoothInterpolation(targetLevel, fine);
}

void MxGeoMultigridPrec::coarsen(int targetLevel, const Epetra_MultiVector & fine, Epetra_MultiVector & coarse) const {
  // smooth before coarsening? (this is the way ML does it...)
  Epetra_MultiVector fineCopy(fine);
  //smoothInterpolation(targetLevel - 1, fineCopy);

  rhsCoarseners[targetLevel - 1]->Apply(fineCopy, coarse);
  if (rcf) MxUtil::Epetra::removeConstField(coarse);
}

double MxGeoMultigridPrec::getResidual(const Epetra_CrsMatrix & op, const Epetra_MultiVector & x, 
const Epetra_MultiVector & b) const {
  Epetra_MultiVector work(x);
  op.Apply(x, work);
  work.Update(-1., b, 1.);
  double res;
  work.Norm2(&res);
  return res;
}

#if 0
void PhcGeoMGPrec::GetBoundaryBulkResiduals (const Epetra_CrsMatrix & op, 
                                             const Epetra_MultiVector & x, 
                                             const Epetra_MultiVector & b,
                                             double & bndryRes,
                                             double & bulkRes,
                                             double & totalRes,
                                             int level) const {
  Epetra_MultiVector work(x);
  op.Apply(x, work);
  work.Update(-1., b, 1.);
  bndryRes=0.; bulkRes=0.; totalRes=0.;

  int numInds = x.Map().getNodeNumElements();
  double * vals = work[0];
  for (int i=0; i<numInds; i++) {
    double val = vals[i];
    double afrac = (*dmAreas_[level])[i];
    if ((afrac==0.) || (afrac==1.))
      bulkRes += val*val;
    else {
      //std::cout << "afrac: " << afrac << std::endl;
      bndryRes += val*val;
    }
    totalRes += val*val;
  }
  bulkRes = sqrt(bulkRes);
  bndryRes = sqrt(bndryRes);
  totalRes = sqrt(totalRes);
}
#endif

int MxGeoMultigridPrec::ApplyInverse(const Epetra_MultiVector & b, Epetra_MultiVector & x) const {
  
  numApplyInverse++;
  time_t start, end;
  time(&start);

  // Allocate workspace for multigrid cycles
  int nvec = b.NumVectors();
  std::vector<Epetra_MultiVector *> bvecs, xvecs;

  //std::cout << *(ops_[levels_-1]);
  //std::cout << x;

  // x will be modified, so put it in workspace
  xvecs.push_back(&x);
  //xvecs.push_back(new Epetra_MultiVector(x));
  // always copy b
  bvecs.push_back(new Epetra_MultiVector(b));
  //RemoveConstField(*bvecs[0]);

  // now allocate the coarse grid stuff, coarsening b all the way down to 
  // initialize the start of a full multigrid cycle
  for (int i = 1; i < levels; ++i) {
    bvecs.push_back(new Epetra_MultiVector(ops[i]->RangeMap(), nvec));
    //coarseners[i - 1]->Apply(*bvecs[i - 1], *bvecs[i]); //coarsen b
    coarsen(i, *bvecs[i - 1], *bvecs[i]);
    //SmoothInterpolation(i, *bvecs[i]);
    xvecs.push_back(new Epetra_MultiVector(ops[i]->DomainMap(), nvec, true)); //zero-out
  }

  // do it!
  fullVCycle(bvecs, xvecs);

  // delete coarsegrids' workspace
  for (int i = 1; i < levels; ++i) {
    delete bvecs[i];
    delete xvecs[i];
  }
  // delete copied b
  delete bvecs[0];
  //x = *xvecs[0];
  //delete xvecs[0];

  //throw 1;
  
  time(&end);
  totalTime += difftime(end, start);

  return 0;
}

int MxGeoMultigridPrec::fullVCycle(std::vector<Epetra_MultiVector*> & bvecs, 
std::vector<Epetra_MultiVector*> & xvecs) const {

  // vars for output
  //double bndryRes, bulkRes, totalRes;
  //char name[200];

  // get coarsest multivectors
  Epetra_MultiVector * x = xvecs[levels - 1];
  Epetra_MultiVector * b = bvecs[levels - 1];

#if 0
  if (output_) {
    spaces(levels_-1); std::cout << "before coarse solve:\n";
    GetBoundaryBulkResiduals(*(ops_[levels_-1]), *x, *b, bndryRes, bulkRes, totalRes, levels_-1);
    spaces(levels_-1); std::cout << "  total residual: " << totalRes << "\n";
    spaces(levels_-1); std::cout << "  bndry residual: " << bndryRes << "\n";
    spaces(levels_-1); std::cout << "  bulk  residual: " << bulkRes << "\n";
  }
#endif

  // solve on coarsest
  //spaces(levels - 1); std::cout << "residual before smooth: " << getResidual(*ops[levels - 1], *x, *b) << "\n";
  smoothers[levels - 1]->ApplyInverse(*b, *x);
  //spaces(levels - 1); std::cout << "residual after smooth: " << getResidual(*ops[levels - 1], *x, *b) << "\n";

#if 0
  if (output_) {
    spaces(levels_-1); std::cout << "after coarse solve:\n";
    GetBoundaryBulkResiduals(*(ops_[levels_-1]), *x, *b, bndryRes, bulkRes, totalRes, levels_-1);
    spaces(levels_-1); std::cout << "  total residual: " << totalRes << "\n";
    spaces(levels_-1); std::cout << "  bndry residual: " << bndryRes << "\n";
    spaces(levels_-1); std::cout << "  bulk  residual: " << bulkRes << "\n";

    spaces(levels_-1); std::cout << "saving coarsest fields\n";
    sprintf(name, "geoMG-fmv-coarse-x.h5");
    bfields_[levels_-1]->SaveToH5(name, comm_, x);
    sprintf(name, "geoMG-fmv-coarse-b.h5");
    bfields_[levels_-1]->SaveToH5(name, comm_, b);
    Epetra_MultiVector dmAinvb(*x);
    dmAinvb.ReciprocalMultiply(1., *dmAreas_[levels_-1], *b, 0.);
    bfields_[levels_-1]->SaveToH5("geoMG-fmv-coarse-dmAinvb.h5", comm_, &dmAinvb);
  }
#endif

  // do walk-up
  for (int i = levels - 2; i >= 0; --i) {
    if (output) {
      spaces(i);
      std::cout << "Interpolating to next FMV level\n";
    }

    //interpolate to current level 
    //xvecs[i+1]->ReciprocalMultiply(1., *dmAreas_[i+1], *xvecs[i+1], 0.);
    //refiners[i]->Apply(*xvecs[i + 1], *xvecs[i]);
    refine(i, *xvecs[i + 1], *xvecs[i]);
    //SmoothInterpolation(i, *xvecs[i]);
    //xvecs[i]->Multiply(1., *dmAreas_[i], *xvecs[i], 0.);

    // Start regular VCycle from current level
    for (int j = 0; j < cycles; ++j) {
      if (output) {
        spaces(i);
        std::cout << "Starting VCycle " << j << ":\n";
      }
      vCycle(bvecs, xvecs, i, j);
    }
  }
  return 0;
}

// for printing
void MxGeoMultigridPrec::spaces(int level) const {
  for (int ii = level; ii < levels - 1; ++ii) std::cout << "  ";
  std::cout << "Level " << level << ": ";
}
