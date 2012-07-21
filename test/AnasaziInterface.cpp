#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include <complex>

#include "MxTypes.h"
#include "MxComm.hpp"
#include "MxMultiVector.hpp"
#include "MxCrsMatrix.hpp"
#include "MxAnasaziMV.hpp"

#include "AnasaziMultiVec.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziMVOPTester.hpp"
#include "AnasaziBasicOutputManager.hpp"

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
#endif

  RCP<MxComm> comm = rcp(new MxComm());

  typedef std::complex<double> Scalar;
  //typedef double Scalar;

  typedef Anasazi::MultiVec<Scalar> MV;
  typedef Anasazi::MultiVecTraits<Scalar, MV> MVT;

  RCP<Anasazi::OutputManager<Scalar> > om =
      rcp(new Anasazi::BasicOutputManager<Scalar>());

  om->setVerbosity(Anasazi::Warnings);


  // create a map
  size_t len = 5;
  size_t nvecs = 2;
  
  RCP<MxMap> map = rcp(new MxMap(len, comm));

  RCP<MxAnasaziMV<Scalar> > ivec =
      rcp(new MxAnasaziMV<Scalar>(map, nvecs));

  bool ierr = Anasazi::TestMultiVecTraits<Scalar, MV>(om, ivec);
  if (ierr) {
    om->print(Anasazi::Warnings,
        "*** MyMultiVec<complex> PASSED TestMultiVecTraits()\n");
  }
  else {
    om->print(Anasazi::Warnings,
        "*** MyMultiVec<complex> FAILED TestMultiVecTraits() ***\n\n");
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
