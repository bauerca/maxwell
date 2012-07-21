#include "MxIO.h"

#include "Epetra_Import.h"

template<size_t DIM, typename Scalar>
MxIO<DIM, Scalar>::MxIO(RCP<MxComm> theComm) : mComm(theComm) {}

template<size_t DIM, typename Scalar>
void MxIO<DIM, Scalar>::openFile(std::string name, char rw) {
  plist_id = H5Pcreate(H5P_FILE_ACCESS);

#ifdef HAVE_MPI
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif

  int pid = mComm->myPID();
  if (pid == 0) std::cout << "Writing " << name << "...\n";

  if (rw == 'w')
    file_id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  else if (rw == 'r')
    file_id = H5Fopen(name.c_str(), H5F_ACC_RDONLY, plist_id);

  H5Pclose(plist_id);
}

template<size_t DIM, typename Scalar>
void MxIO<DIM, Scalar>::closeFile() {
  H5Fclose(file_id);
}

template<size_t DIM, typename Scalar>
void MxIO<DIM, Scalar>::saveLinear(int numComps, int globLen, int myLen, int myOffset, double * data, 
std::string name) {

  openFile(name, 'w');
  
  int dset_dims = (numComps > 1 ? 2 : 1);

  // create the dataspace for the dataset
  dimsf[0] = globLen;
  if (dset_dims == 2) dimsf[1] = numComps;
  filespace = H5Screate_simple(dset_dims, dimsf, NULL);

  // for hdf5 >= 1.8
#if H5_VERS_MINOR == 8
  llist_id = H5Pcreate(H5P_LINK_CREATE);
  alist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dcreate2(file_id, "data", H5T_NATIVE_DOUBLE, filespace, llist_id, H5P_DEFAULT, alist_id);
  H5Pclose(llist_id);
  H5Pclose(alist_id);
#else
  // for hdf5 1.6
  dset_id = H5Dcreate(file_id, "data", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT);
#endif

  //H5Sclose(filespace);

  chunk_dims[0] = myLen;
  count[0] = 1;
  offset[0] = myOffset;
  stride[0] = 1;
  block[0] = chunk_dims[0];
  if (dset_dims == 2) {
    chunk_dims[1] = numComps;
    count[1] = 1;
    offset[1] = 0;
    stride[1] = 1;
    block[1] = chunk_dims[1];
  }

  memspace  = H5Screate_simple(dset_dims, chunk_dims, NULL);

  // Select hyperslab in the file
  //filespace = H5Dget_space(dset_id);
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif

  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);

  // Close/release resources
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(plist_id);
  H5Fclose(file_id);
}

template<size_t DIM, typename Scalar>
void MxIO<DIM, Scalar>::saveOnDomain(int numComps,
MxGrid<DIM> const & theGrid, double * data, std::string filename,
std::string datasetName) {
  //openFile(name, 'w');

  int dset_dim = DIM + 1;
  // create the dataspace for the dataset
  MxGridDomain<DIM> domain = theGrid.getGridDomain(0); // zero guard cells
  MxDimVector<int, DIM> gridRes = theGrid.getResolution() + MxDimVector<int, DIM>(1);
  MxDimVector<int, DIM> domainRes = domain.getInteriorResolution();

  for (size_t i = 0; i < dset_dim; ++i) {
    if (i == DIM) {
      dimsf[i] = numComps;
      chunk_dims[i] = numComps;
    }
    else {
      dimsf[i] = gridRes[i];
      chunk_dims[i] = domainRes[i];
    }
  }

  MxDimVector<int, DIM> lowerBound = domain.getLowerBoundCell();

  for (size_t i = 0; i < dset_dim; ++i) {
    if (i == DIM)
      offset[i] = 0;
    else
      offset[i] = lowerBound[i];
    count[i] = 1;
    stride[i] = 1;
    block[i] = chunk_dims[i];
  }

  filespace = H5Screate_simple(dset_dim, dimsf, NULL);
  memspace  = H5Screate_simple(dset_dim, chunk_dims, NULL);

  // for hdf5 >= 1.8
#if H5_VERS_MINOR == 8
  llist_id = H5Pcreate(H5P_LINK_CREATE);
  alist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dcreate2(file_id, datasetName.c_str(), H5T_NATIVE_DOUBLE,
      filespace, llist_id, H5P_DEFAULT, alist_id);
  H5Pclose(llist_id);
  H5Pclose(alist_id);
#else
  // for hdf5 1.6
  dset_id = H5Dcreate(file_id, datasetName.c_str(), H5T_NATIVE_DOUBLE,
      filespace, H5P_DEFAULT);
#endif

  H5Sclose(filespace);

  // Select hyperslab in the file
  filespace = H5Dget_space(dset_id);
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);
  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif

  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);

  // Close/release resources
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(plist_id);
  //H5Fclose(file_id);
}

template<size_t DIM, typename Scalar>
void MxIO<DIM, Scalar>::save(MxMultiVector<Scalar> const & fieldData,
MxGridField<DIM> const & theGridField, std::string name) {
  if (*fieldData.getMap() != *theGridField.getMap()) {
    std::cout << "MxIO<DIM, Scalar>::save(fieldData, gridField, name): dist map of fieldData does not match that of gridField. Field save aborted.\n";
    return;
  }

  int numComps = theGridField.getNumComps();

  MxGrid<DIM> const & grid = theGridField.getGrid();
  MxGridDomain<DIM> domain = grid.getGridDomain(0);
  int N = domain.getNumInteriorCells();
  int numVals = numComps * N;
  Scalar * data = new Scalar[numComps * N];

  char filename[200];
  std::string fieldName(theGridField.getName());
  MxDimVector<int, DIM> cell;
  MxIndex lid;
  for (int vec = 0; vec < fieldData.getNumVecs(); ++vec) {
    sprintf(filename, "mx%s_%s_vec%.2i.h5", fieldName.c_str(), name.c_str(), vec);

    int j = 0;
    for (int i = 0; i < N; ++i) {
      cell = domain.interiorIndxToCell(i);
      for (size_t comp = 0; comp < numComps; ++comp) {
        lid = fieldData.getMap()->getLocalIndex(
            theGridField.globCompIndx(comp, cell));
        //std::cout << lid << "\n";
        if (lid != MxInvalidIndex)
          data[j] = fieldData(lid, vec);
        else
          data[j] = 0;
        j++;
      }
    }

    openFile(filename, 'w');
    if (ScalarTraits<Scalar>::isComplex) {
      double * datapart = new double[numVals];
      for (size_t i = 0; i < numVals; ++i)
        datapart[i] = double(ScalarTraits<Scalar>::real(data[i]));
      saveOnDomain(numComps, grid, datapart, filename, std::string("real"));
      for (size_t i = 0; i < numVals; ++i)
        datapart[i] = double(ScalarTraits<Scalar>::imag(data[i]));
      saveOnDomain(numComps, grid, datapart, filename, std::string("imag"));
      delete[] datapart;
    }
    else {
      saveOnDomain(numComps, grid, (double *) data, filename, std::string("real"));
    }
    closeFile();
  }

  delete[] data;
}

#if 0
template<size_t DIM, typename Scalar>
void MxIO<DIM, Scalar>::save(std::vector<MxDimVector<double, DIM> > const & points,
MxMultiVector<Scalar> const & fieldData, MxGridField<DIM> const & theGridField, std::string name) {
  
  MxMultiVector<Scalar> * fieldValues;
  int numComps = theGridField.getNumComps();

  char filename[200];
  std::string fieldName(theGridField.getName());
  fieldValues = theGridField.getFieldValues(points, fieldData);
  for (int vec = 0; vec < fieldData.NumVectors(); ++vec) {
    sprintf(filename, "mx%s_%s_vec%.2i.h5", fieldName.c_str(), name.c_str(), vec);
      
    saveLinear(numComps, fieldValues->GlobalLength() / numComps, fieldValues->MyLength()/ numComps, fieldValues->Map().GID(0) - fieldValues->Map().IndexBase(), (*fieldValues)[vec], filename);
  }
  delete fieldValues;
}

template<size_t DIM, typename Scalar>
void MxIO<DIM, Scalar>::save(MxPointCloud<DIM> const & pointCloud,
MxMultiVector<Scalar> const & fieldData, MxGridField<DIM> const & theGridField, std::string name) {
  
  MxMultiVector<Scalar> * fieldValues, * linFieldVals;
  int numComps = theGridField.getNumComps();

  char filename[200];
  std::string fieldName(theGridField.getName());
  fieldValues = theGridField.getFieldValues(pointCloud, fieldData);

  // import values if fieldValues is not on a linear map
  if (fieldValues->Map().LinearMap())
    linFieldVals = fieldValues;
  else {
    MxMap linearMap(*pointCloud.getLinearFieldMap(numComps));
    linFieldVals = new MxMultiVector<Scalar>(linearMap, fieldValues->NumVectors());
    Epetra_Import importer(linearMap, fieldValues->Map());
    linFieldVals->Import(*fieldValues, importer, Insert);
  }
  int gLen = linFieldVals->GlobalLength() / numComps;
  int mLen = linFieldVals->MyLength() / numComps;
  int offs = (linFieldVals->Map().GID(0) - linFieldVals->Map().IndexBase()) / numComps;

  // save the data
  for (int vec = 0; vec < fieldData.NumVectors(); ++vec) {
    sprintf(filename, "mx%s_%s_vec%.2i.h5", fieldName.c_str(), name.c_str(), vec);
    saveLinear(numComps, gLen, mLen, offs, (*linFieldVals)[vec], filename);
  }

  // delete data
  if (!fieldValues->Map().LinearMap())
    delete linFieldVals;
  delete fieldValues;
}

template<size_t DIM, typename Scalar>
void MxIO<DIM, Scalar>::save(MxGrid<DIM> const & theGrid, MxMultiVector<Scalar> const & fieldData,
MxGridField<DIM> const & theGridField, std::string name) {
  if (theGrid == theGridField.getGrid()) {
    save(fieldData, theGridField, name);
  }
  else {
    std::cout << "Saving field values on an arbitrary grid not implemented yet. Soon!\n";
  }
}

template<size_t DIM, typename Scalar>
void MxIO<DIM, Scalar>::load(std::string filename, MxGridField<DIM> * theGridField) {

  openFile(filename, 'r');


}

#endif



template class MxIO<1, double>;
template class MxIO<2, double>;
template class MxIO<3, double>;
template class MxIO<1, MxComplex>;
template class MxIO<2, MxComplex>;
template class MxIO<3, MxComplex>;
