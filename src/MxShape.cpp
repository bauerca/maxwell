
#include "MxShape.hpp"
#include "MxUtil.hpp"
#include "MxGrid.h"
#include "MxGridDomainIter.hpp"

#include "hdf5.h"

#include "Epetra_Comm.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#endif

template<size_t DIM>
MxDimVector<double, DIM> MxShape<DIM>::boundaryPoint(const MxDimVector<double, DIM> & p1, const MxDimVector<double, DIM> & p2) const {
  return MxUtil::rootFind<MxShape<DIM>, DIM>(*this, p1, p2, 1.0e-12, 20);
}

template<size_t DIM>
void MxShape<DIM>::applyTransforms(MxShape<DIM> & shape) const {
  shape.A = A * shape.A;
  shape.Ainv = shape.A.inv();
  shape.sign = sign;
  shape.b += b;
  shape.Ainvb = shape.Ainv * shape.b;
}

template<size_t DIM>
void MxShape<DIM>::transforms(const Teuchos::XMLObject & node) {
  // should we turn shape inside out?
  std::string io = MxUtil::XML::getAttr("invert", node, "false");
  if (io == "true")
    sign = -1.0;
  else
    sign = 1.0;

  int n = node.numChildren();
  Teuchos::XMLObject child;
  std::string tag;
  for (int i = 0; i < n; ++i) {
    child = node.getChild(i);
    tag = child.getTag();

    if (tag == "Rotate")
      rotate(child);
    else if (tag == "Translate")
      translate(child);
    else if (tag == "Reflect")
      reflect(child);
    else if (tag == "Scale")
      scale(child);
    else {
      //std::cout << "MxShape::transforms(...): Warning: did not understand tag, '"
                //<< tag << "' in XML block:\n";
      //child.print(std::cout, 1);
      continue;
    }
  }
}

template<size_t DIM>
void MxShape<DIM>::rotate(const Teuchos::XMLObject & node) {
  MxDimVector<double, 3> pivot, axis;

  axis.strFill(MxUtil::XML::getAttr("axis", node));
  double angle = atof(MxUtil::XML::getAttr("angle", node).c_str());
  std::string s = MxUtil::XML::getAttr("pivot", node, "");
  if (s != "") {
    pivot.strFill(s);
    rotate(axis, angle, pivot);
  }
  else
    rotate(axis, angle);
}

/*
template<size_t DIM>
void MxShape<DIM>::rotate(const MxDimVector<double, 1> & axis, double angle) {
  MxDimVector<double, 3> a(axis / axis.norm());
  rotate(a, angle);
}

template<size_t DIM>
void MxShape<DIM>::rotate(const MxDimVector<double, 2> & axis, double angle) {
  MxDimVector<double, 3> a(axis / axis.norm());
  rotate(a, angle);
}
*/

template<size_t DIM>
void MxShape<DIM>::rotate(const MxDimVector<double, 3> & axis, double angle) {
  // this rotates about the current translation point
  rotate(axis, angle, b);
  //rotate(axis, angle, 0);
}

template<size_t DIM>
void MxShape<DIM>::rotate(const MxDimVector<double, 3> & axis, double angle, const MxDimVector<double, 3> & pivot) {
  // normalize 'axis'. Conversion from smaller DimVecs to DimVec3 preserves 
  // unity because the default conversion fills the vector with zeros.
  MxDimVector<double, 3> a(axis / axis.norm());

  // R is the rotation matrix. M is the cross-product matrix for 'axis'
  MxDimMatrix<double, 3> R, M, eye(MxDimVector<double, 3>(1.0));

  // fill M
  M(0, 1) = M(1, 0) = a[2];  M(1, 0) *= -1.0; 
  M(2, 0) = M(0, 2) = a[1];  M(0, 2) *= -1.0;
  M(1, 2) = M(2, 1) = a[0];  M(2, 1) *= -1.0;

  // find R using Rodrigues' rotation formula
  R = eye + M * sin(angle) + (M * M) * (1.0 - cos(angle));

  // adjust transformation
  A = R * A;
  b = R * (b - pivot) + pivot;

  //std::cout << "Linear transformation: A\n";
  //A.print();
  //std::cout << "Translation: b\n";
  //b.print();

  Ainv = A.inv();
  Ainvb = Ainv * b;
}

template<size_t DIM>
void MxShape<DIM>::scale(const Teuchos::XMLObject & node) {
  MxDimVector<double, 3> magnitudes;
  // if user doesn't provide a vector of 3 magnitudes, fill the rest with
  // with 1s
  magnitudes.strFill(MxUtil::XML::getAttr("magnitudes", node), 1.0); 
  scale(magnitudes);
}

/*
template<size_t DIM>
void MxShape<DIM>::scale(const MxDimVector<double, 1> & magnitudes) {
  MxDimVector<double, 3> mags3(magnitudes, 1.0); // pad with 1s
  scale(mags3);
}

template<size_t DIM>
void MxShape<DIM>::scale(const MxDimVector<double, 2> & magnitudes) {
  MxDimVector<double, 3> mags3(magnitudes, 1.0); // pad with 1s
  scale(mags3);
}
*/

template<size_t DIM>
void MxShape<DIM>::scale(const MxVecD3 & magnitudes) {
  scale(magnitudes, 0);
}

template<size_t DIM>
void MxShape<DIM>::scale(const MxDimVector<double, 3> & magnitudes, const MxDimVector<double, 3> & origin) {
  MxDimMatrix<double, 3> S(magnitudes);

  // adjust transformation
  A = S * A;
  b = S * (b - origin) + origin;
  Ainv = A.inv();
  Ainvb = Ainv * b;
}


template<size_t DIM>
void MxShape<DIM>::translate(const Teuchos::XMLObject & node) {
  MxDimVector<double, 3> v;
  v.strFill(MxUtil::XML::getAttr("vector", node));
  translate(v);
}

template<size_t DIM>
void MxShape<DIM>::translate(const MxVecD3 & v) {
  // adjust transformation
  b += v;

  Ainvb = Ainv * b;
}

template<size_t DIM>
void MxShape<DIM>::reflect(const Teuchos::XMLObject & node) {
  MxDimVector<double, 3> n, p;
  n.strFill(MxUtil::XML::getAttr("normal", node));
  p.strFill(MxUtil::XML::getAttr("point in plane", node));
  reflect(n, p);
}

template<size_t DIM>
void MxShape<DIM>::reflect(const MxDimVector<double, 3> & normal, const MxDimVector<double, 3> & pointInPlane) {
  // normalize normal
  MxDimVector<double, 3> n(normal / normal.norm());

  // M is the reflection about the plane that intersects the origin.
  MxDimMatrix<double, 3> N(n, n), eye(MxDimVector<double, 3>(1.0)), M;
  M = eye - 2.0 * N;

  // adjust transformation
  A = M * A;
  b = M * (b - pointInPlane) + pointInPlane;

  Ainv = A.inv();
  Ainvb = Ainv * b;
}
  

template<size_t DIM>
void MxShape<DIM>::save(const MxGrid<DIM> & aGrid, const char * dSetName, const char * fileName) const {

  /*
   * HDF5 APIs definitions
  */   

  hid_t       file_id, dset_id;         /* file and dataset identifiers */
  hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
  hsize_t     dimsf[DIM];                 /* dataset dimensions */
  hsize_t     chunk_dims[DIM];            /* chunk dimensions */
  hsize_t count[DIM];           /* hyperslab selection parameters */
  hsize_t stride[DIM];
  hsize_t block[DIM];
  hsize_t offset[DIM];
  hid_t plist_id;                 /* property list identifier */
  hid_t llist_id;                 /* link creation list identifier */
  hid_t alist_id;                 /* dataset access list identifier */
  herr_t  status;

  char filename[200];
  if (fileName != 0)
    sprintf(filename, "mxShapeFunc_%s.h5", fileName);
  else
    sprintf(filename, "mxShapeFunc_%s.h5", &name[0]);

  plist_id = H5Pcreate(H5P_FILE_ACCESS);

#ifdef HAVE_MPI
  //MPI_Comm * comm = &dynamic_cast<const Epetra_MpiComm &>(aGrid.getComm()).Comm();
  //MPI_Comm comm = MPI_COMM_WORLD;
  //MPI_Info info = MPI_INFO_NULL;

  //H5Pset_fapl_mpio(plist_id, comm, info);
  //
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif

  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  // create the dataspace for the dataset
  MxGridDomain<DIM> domain = aGrid.getGridDomain(0); // zero guard cells
  MxDimVector<int, DIM> gridRes = aGrid.getResolution() + MxDimVector<int, DIM>(1);
  MxDimVector<int, DIM> domainRes = domain.getInteriorResolution();

  for (size_t i = 0; i < DIM; ++i) {
    dimsf[i] = gridRes[i];
    chunk_dims[i] = domainRes[i];
  }

  filespace = H5Screate_simple(DIM, dimsf, NULL);
  memspace  = H5Screate_simple(DIM, chunk_dims, NULL);

  // create chunked dataset
  //plist_id = H5Pcreate(H5P_DATASET_CREATE);
  //H5Pset_chunk(plist_id, DIM, chunk_dims);

  // for hdf5 >= 1.8
#if H5_VERS_MINOR == 8
  llist_id = H5Pcreate(H5P_LINK_CREATE);
  alist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dcreate2(file_id, dSetName, H5T_NATIVE_DOUBLE, filespace, llist_id, H5P_DEFAULT, alist_id);
  H5Pclose(llist_id);
  H5Pclose(alist_id);
#else
  // for hdf5 1.6
  dset_id = H5Dcreate(file_id, dSetName, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT);
#endif

  H5Sclose(filespace);

  MxDimVector<int, DIM> lowerBound = domain.getLowerBoundCell();
  
  for (size_t i = 0; i < DIM; ++i) {
    count[i] = 1;
    offset[i] = lowerBound[i];
    stride[i] = 1;
    block[i] = chunk_dims[i];
  }

  // Select hyperslab in the file
  filespace = H5Dget_space(dset_id);
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);

  // set data
  int N = domain.getNumInteriorCells();
  double * data = new double[N];
  MxDimVector<double, DIM> node;
  MxDimVector<int, DIM> cell;
  for (int i = 0; i < N; ++i) {
    cell = domain.interiorIndxToCell(i);
    node = aGrid.nodeCoord(cell);
    data[i] = this->func(node);
  }

  // Create property list for collective dataset write.
  plist_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif

  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);

  delete[] data;

  // Close/release resources
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(plist_id);
  H5Fclose(file_id);

}


template class MxShape<1>;
template class MxShape<2>;
template class MxShape<3>;
