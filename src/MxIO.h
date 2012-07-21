#ifndef MX_IO
#define MX_IO

#include <string>
#include <vector>

#include "MxDimVector.hpp"
#include "MxGrid.h"
#include "MxGridField.hpp"
#include "MxShape.hpp"
#include "MxPointCloud.h"

#include "hdf5.h"

#include "MxComm.hpp"
#include "MxMultiVector.hpp"

template<size_t DIM, typename Scalar>
class MxIO {
  public:
    MxIO(RCP<MxComm> comm);

    void save(MxMultiVector<Scalar> const & fieldData,
        MxGridField<DIM> const & theGridField, std::string name);

    void save(std::vector<MxDimVector<double, DIM> > const & points,
        MxMultiVector<Scalar> const & fieldData,
        MxGridField<DIM> const & theGridField, std::string name);

    void save(MxPointCloud<DIM> const & pointCloud,
        MxMultiVector<Scalar> const & fieldData,
        MxGridField<DIM> const & theGridField, std::string name);

    void save(MxGrid<DIM> const & theGrid,
        MxMultiVector<Scalar> const & fieldData,
        MxGridField<DIM> const & theGridField, std::string name);

    void save(MxShape<DIM> const & theShape, MxGrid<DIM> const & theGrid);
    
    
    void load(std::string filename,
        std::vector<MxDimVector<double, DIM> > & points);

    void load(std::string filename, MxGridField<DIM> * theGridField);

  private:

    void openFile(std::string name, char rw);

    void closeFile();

    void saveOnDomain(int numComps, MxGrid<DIM> const & theGrid,
        double * data, std::string filename, std::string datasetName);

    void saveLinear(int numComps, int globLen, int myLen, int myOffset,
        double * data, std::string name);
    
    RCP<MxComm> mComm;

  /*
   * HDF5 APIs definitions
  */
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dimsf[4];                 /* dataset dimensions */
    hsize_t     chunk_dims[4];            /* chunk dimensions */
    hsize_t count[4];           /* hyperslab selection parameters */
    hsize_t stride[4];
    hsize_t block[4];
    hsize_t offset[4];
    hid_t plist_id;                 /* property list identifier */
    hid_t llist_id;                 /* link creation list identifier */
    hid_t alist_id;                 /* dataset access list identifier */
    herr_t  status;

};

#endif
