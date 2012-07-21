#ifndef MX_INTERPOLATOR
#define MX_INTERPOLATOR

#include <vector>

#include "MxTypes.h"
#include "MxDimVector.hpp"
#include "MxGridField.hpp"
#include "MxPointCloud.h"
#include "MxMap.hpp"
#include "MxCrsMatrix.hpp"

template<size_t DIM, typename Scalar>
class MxGridFieldInterpolator : public MxCrsMatrix<Scalar> {
  public:
    MxGridFieldInterpolator(MxGridField<DIM> const * fromField,
      MxGridField<DIM> const * toField);

    MxGridFieldInterpolator(MxGridField<DIM> const * fromField,
      RCP<MxPointCloud<DIM> const> pointCloud);

    virtual ~MxGridFieldInterpolator() {};

    //virtual const MxGridField<DIM> * getDomainField() const {return field;}

    //virtual const MxGridField<DIM> * getRangeField() const {return field;}

    //virtual const MxGrid<DIM> * getGrid() const {return grid;}

  protected:

    virtual void getStencil(size_t comp, MxDimVector<double, DIM> point,
      std::vector<MxIndex> & indices, std::vector<Scalar> & vals);

  private:

    MxGridField<DIM> const * mFromField, * mToField;

    RCP<MxMap> mFromMap, mToMap;

    RCP<MxPointCloud<DIM> const> mPointCloud;

    void setMatrixField();
    
    void setMatrixCloud();

    void insertValues(std::vector<MxIndex> const & rows,
      std::vector<std::vector<MxIndex> > const & cols,
      std::vector<std::vector<Scalar> > const & vals);

};

#endif
