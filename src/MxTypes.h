#ifndef MX_TYPES
#define MX_TYPES

#include <complex>

#include "Teuchos_RCP.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ArrayView.hpp"

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#define EPETRA 0
#define TPETRA 1

#define LINALG_BASE EPETRA

enum MxPolType {TM, TE};

enum MxBCType {
  ZERO, // Dirichlet with value set to zero
  CONSTANT, // Neumann with normal derivative set to zero
  PERIODIC,
  PEC,
  PMC
};

// Scalar type enumeration
enum MxScalarType {
  Complex,
  Real
};

namespace Mx {
  using Teuchos::RCP;
  using namespace Tpetra;
}

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::rcp;
using Teuchos::ScalarTraits;


// Complex type for everything
typedef std::complex<double> MxComplex;

// Ordinal type for everything
typedef size_t MxIndex;
const MxIndex MxInvalidIndex = Teuchos::OrdinalTraits<MxIndex>::invalid();


// annoying views of arrays required by Tpetra
typedef Teuchos::ArrayView<MxIndex> MxIndexView;
typedef Teuchos::ArrayView<size_t> MxSizeTView;
typedef Teuchos::ArrayView<int> MxIntView;


#endif
