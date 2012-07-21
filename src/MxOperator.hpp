#ifndef MX_OPERATOR
#define MX_OPERATOR

//template<size_t DIM> class MxGridField;
//template<size_t DIM> class MxGrid;
//

#include "MxGridField.hpp"
#include "MxGrid.h"

template<size_t DIM>
class MxOperator {
  public:
    //virtual void setMatrix() = 0;

    //virtual const Epetra_CrsMatrix * getMatrixPtr() const = 0;

    //virtual Epetra_CrsMatrix getMatrixCopy() const = 0;

    //virtual void clearMatrix() = 0;

    virtual ~MxOperator() {}

    virtual const MxGridField<DIM> * getDomainField() const = 0;

    virtual const MxGridField<DIM> * getRangeField() const = 0;

    virtual const MxGrid<DIM> * getGrid() const = 0;

    virtual void setImaginary(bool imag) {imaginary = imag;}

    virtual bool isImaginary() {return imaginary;}

  protected:
    bool imaginary;

};


#endif
