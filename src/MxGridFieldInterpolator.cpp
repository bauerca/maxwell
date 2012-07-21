#include "MxGridFieldInterpolator.h"
#include "MxDimVector.hpp"
#include "MxGridFieldIter.hpp"


template<size_t DIM, typename Scalar>
MxGridFieldInterpolator<DIM, Scalar>::MxGridFieldInterpolator(
MxGridField<DIM> const * fromField, 
MxGridField<DIM> const * toField) :
MxCrsMatrix<Scalar>(toField->getMap(), fromField->getMap()),
mToField(toField), mFromField(fromField),
mToMap(toField->getMap()), mFromMap(fromField->getMap()) {
  setMatrixField();
}

template<size_t DIM, typename Scalar>
MxGridFieldInterpolator<DIM, Scalar>::MxGridFieldInterpolator(
MxGridField<DIM> const * fromField, 
RCP<MxPointCloud<DIM> const> pointCloud) :
MxCrsMatrix<Scalar>(pointCloud->getFieldMap(fromField->getNumComps()), fromField->getMap()),
mFromField(fromField), mPointCloud(pointCloud),
mToMap(pointCloud->getFieldMap(fromField->getNumComps())),
mFromMap(fromField->getMap()) {
  setMatrixCloud();
}

template<size_t DIM, typename Scalar>
void MxGridFieldInterpolator<DIM, Scalar>::getStencil(size_t comp,
MxDimVector<double, DIM> point, std::vector<MxIndex> & cols,
std::vector<Scalar> & vals) {
  cols.clear();
  vals.clear();

  MxGrid<DIM> const & fromGrid = mFromField->getGrid();  

  MxDimVector<double, DIM> cellCoord(mFromField->getCompCellCoord(comp));
  MxDimVector<double, DIM> origin(fromGrid.getOrigin());
  MxDimVector<double, DIM> cellSize(fromGrid.getCellSize());
  MxDimVector<double, DIM> p, p0, ones(1);

  MxIndex col;
  Scalar val;

  // find the nearest component of 'this' field whose coordinates are
  // all lower than those of p.
  MxDimVector<int, DIM> myCell;
  for (size_t i = 0; i < DIM; i++)
    myCell[i] = int(floor((point[i] - cellCoord[i] - origin[i]) /
        cellSize[i]));

  std::vector<MxDimVector<int, DIM> > cellBlock;
  cellBlock = fromGrid.getCellBlock(myCell, 2);

  for (size_t i = 0; i < cellBlock.size(); i++) {
    //cellBlock[i].print();
    p0 = fromGrid.nodeCoord(cellBlock[i]) + cellCoord;
    col = mFromField->globCompIndx(comp, cellBlock[i]);

    // implicit cast to complex if required
    val = (ones - fabs(p - p0) / cellSize).prod();

    cols.push_back(col);
    vals.push_back(val);
  }
}

template<size_t DIM, typename Scalar>
void MxGridFieldInterpolator<DIM, Scalar>::insertValues(
std::vector<MxIndex> const & rows,
std::vector<std::vector<MxIndex> > const & indices,
std::vector<std::vector<Scalar> > const & values) {
  size_t numRows = indices.size();
  size_t numEntries;
  MxIndex row, col;
  Scalar val;

  std::vector<MxIndex> offProcRows, offProcCols;
  std::vector<Scalar> offProcVals;

  for (size_t i = 0; i < numRows; ++i) {
    row = rows[i];
    numEntries = indices[i].size(); 
    for (size_t j = 0; j < numEntries; ++j) {
      col = indices[i][j];
      val = values[i][j];
      if (mFromMap->isNodeGlobalIndex(col)) {
        MxCrsMatrix<Scalar>::insertRowValues(row, 1, &col, &val);
      }
      // otherwise if the column GID is off process or outside the
      // field 'region', save the info
      else {
        offProcVals.push_back(val);
        offProcCols.push_back(col);
        offProcRows.push_back(row);
      }
    }
  }

  // find out which global column ids are inside the field 'region'
  // (lid in lidList
  // will be -1 if GID does not exist on any processor)
  std::vector<int> pids(offProcCols.size());
  std::vector<MxIndex> lids(offProcCols.size());
  mFromMap->getRemoteIndexList(offProcCols, pids, lids);

  //MxUtil::printStdVector(cols);

  // fill the off-processor matrix elements
  for (size_t i = 0; i < offProcCols.size(); i++) {
    //std::cout << lidList[i] << ", ";
    if (lids[i] != MxInvalidIndex)
      MxCrsMatrix<Scalar>::insertRowValues(offProcRows[i], 1,
        &offProcCols[i], &offProcVals[i]);
  }
  //std::cout << "\n";

  // scale the rows so that their sums are unity
  MxCrsMatrix<Scalar>::fillComplete(mFromMap, mToMap);

  RCP<MxVector<Scalar> > invRowSums = MxCrsMatrix<Scalar>::getRowSums(true);
  MxCrsMatrix<Scalar>::leftScale(invRowSums);
}

template<size_t DIM, typename Scalar>
void MxGridFieldInterpolator<DIM, Scalar>::setMatrixCloud() {
  size_t numComps = mFromField->getNumComps();

  std::vector<MxIndex> rows(
      mPointCloud->numNodePoints() * numComps);
  std::vector<std::vector<MxIndex> > stencilInds(
      mPointCloud->numNodePoints() * numComps);
  std::vector<std::vector<Scalar> > stencilVals(
      mPointCloud->numNodePoints() * numComps);

  std::vector<MxDimVector<double, DIM> > const & points = 
      mPointCloud->getNodePoints();

  typename std::vector<MxDimVector<double, DIM> >::const_iterator ptIter;
  MxIndex i = 0;
  for (ptIter = points.begin(); ptIter != points.end(); ++ptIter) {
    for (size_t comp = 0; comp < numComps; ++comp) {
      rows[i] = mToMap->getGlobalIndex(i);
      this->getStencil(comp, *ptIter, stencilInds[i], stencilVals[i]);
      i++; 
    }
  }

  insertValues(rows, stencilInds, stencilVals);
}

template<size_t DIM, typename Scalar>
void MxGridFieldInterpolator<DIM, Scalar>::setMatrixField() {
  size_t numRows = mToMap->getNodeNumIndices();
  std::vector<MxIndex> rows(numRows);
  std::vector<std::vector<MxIndex> > stencilInds(numRows);
  std::vector<std::vector<Scalar> > stencilVals(numRows);

  MxGridFieldIter<DIM> iter(mToField);
  size_t i = 0;
  for (iter.begin(); !iter.atEnd(); iter.bump()) {
    rows[i] = iter.getGlobCompIndx();
    this->getStencil(iter.getComp(), iter.getCoord(), stencilInds[i],
        stencilVals[i]);
    i++;
  }

  insertValues(rows, stencilInds, stencilVals);
}


template class MxGridFieldInterpolator<1, double>;
template class MxGridFieldInterpolator<2, double>;
template class MxGridFieldInterpolator<3, double>;
template class MxGridFieldInterpolator<1, MxComplex>;
template class MxGridFieldInterpolator<2, MxComplex>;
template class MxGridFieldInterpolator<3, MxComplex>;
