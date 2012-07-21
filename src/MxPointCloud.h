#ifndef MX_POINT_CLOUD
#define MX_POINT_CLOUD

#include <vector>

#include "MxDimVector.hpp"
#include "MxGrid.h"
#include "MxMap.hpp"

#include "Teuchos_XMLObject.hpp"
#include "Teuchos_RCP.hpp"

template<size_t DIM>
class MxPointCloud {
  public:
    MxPointCloud(std::vector<MxDimVector<double, DIM> > const & points,
      MxGrid<DIM> const * grid);

    MxPointCloud(Teuchos::XMLObject const & node, MxGrid<DIM> const * grid);

    Teuchos::RCP<MxMap> getMap() const {return mMap;}

    Teuchos::RCP<MxMap> getFieldMap(size_t numComps) const;

    Teuchos::RCP<MxMap> getLinearMap() const;

    Teuchos::RCP<MxMap> getLinearFieldMap(int numComps) const;

    std::vector<MxDimVector<double, DIM> > const & getNodePoints() const {return mNodePoints;}

    size_t numNodePoints() const {return mMap->getNodeNumIndices();}

    std::vector<MxIndex> const & getNodePointInds() const {return mNodePtInds;}

    int numPoints() const {return mMap->getGlobalNumIndices();}


  private:
    MxGrid<DIM> const * mGrid;

    std::vector<MxDimVector<double, DIM> > mPoints;

    std::vector<MxDimVector<double, DIM> > mNodePoints;

    std::vector<MxIndex> mNodePtInds;

    Teuchos::RCP<MxMap> mMap;

    void divide();
};

#endif
