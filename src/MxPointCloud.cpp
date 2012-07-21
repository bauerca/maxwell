#include "MxPointCloud.h"

#include "MxUtil.hpp"
#include "MxGridDomain.h"

template<size_t DIM>
MxPointCloud<DIM>::MxPointCloud(
std::vector<MxDimVector<double, DIM> > const & points,
MxGrid<DIM> const * grid) : mGrid(grid), mPoints(points) {
  divide();
}

template<size_t DIM>
MxPointCloud<DIM>::MxPointCloud(Teuchos::XMLObject const & node,
MxGrid<DIM> const * grid) : mGrid(grid) {
  std::string tag = node.getTag();
  if (tag != "Points") {
    std::cout << "MxPointCloud::MxPointCloud(...): cannot create point cloud from a non <Points> block.\n";
    std::cout << "  You have created an empty point cloud.\n";
  }
  else {
    MxDimVector<double, DIM> pt;
    std::vector<std::string> lines(MxUtil::XML::nodeLines(node));
    std::vector<std::string>::const_iterator iter;
    for (iter = lines.begin(); iter != lines.end(); ++iter) {
      pt.strFill(MxUtil::Strings::splitEquation(*iter).second);
      mPoints.push_back(pt);
      //std::cout << *iter << std::endl;
    }
    divide();
  }
}

template<size_t DIM>
RCP<MxMap> MxPointCloud<DIM>::getFieldMap(size_t numComps) const {
  RCP<MxMap> res;

  std::vector<MxIndex> ptFldInds;
  ptFldInds.reserve(numComps * mNodePtInds.size());

  std::vector<MxIndex>::const_iterator iter;
  for (iter = mNodePtInds.begin(); iter != mNodePtInds.end(); ++iter) {
    for (int comp = 0; comp < numComps; ++comp) {
      ptFldInds.push_back(numComps * (*iter) + comp);
    }
  }
  res = rcp(new MxMap(numComps * numPoints(), ptFldInds, mGrid->getComm()));
  return res;
}

template<size_t DIM>
RCP<MxMap> MxPointCloud<DIM>::getLinearMap() const {
  RCP<MxMap> res;
  res = rcp(new MxMap(numPoints(), mGrid->getComm()));
  return res;
}

template<size_t DIM>
RCP<MxMap> MxPointCloud<DIM>::getLinearFieldMap(int numComps) const {
  RCP<MxMap> res;
  
  RCP<MxMap> linPtMap(getLinearMap());
  size_t numMyPts = linPtMap->getNodeNumIndices();

  res = rcp(new MxMap(numComps * numPoints(), numComps * numMyPts, mGrid->getComm()));
  return res;
}

template<size_t DIM>
void MxPointCloud<DIM>::divide() {
  MxGridDomain<DIM> domain(mGrid->getGridDomain(0));

  MxDimVector<double, DIM> lb(mGrid->nodeCoord(domain.getLowerBoundCell()));
  MxDimVector<double, DIM> ub(mGrid->nodeCoord(domain.getUpperBoundCell()));

  typename std::vector<MxDimVector<double, DIM> >::const_iterator iter;

  double xi;
  bool inside;
  for (iter = mPoints.begin(); iter != mPoints.end(); ++iter) {
    inside = true;
    for (size_t i = 0; i < DIM; ++i) {
      xi = (*iter)[i];
      if (xi < lb[i] or xi >= ub[i]) {
        inside = false;
        break;
      }
    }
    if (inside) {
      mNodePtInds.push_back(iter - mPoints.begin());
      mNodePoints.push_back(*iter);
    }
  }

  mMap = rcp(new MxMap(mPoints.size(), mNodePtInds, mGrid->getComm()));
}

template class MxPointCloud<1>;
template class MxPointCloud<2>;
template class MxPointCloud<3>;
