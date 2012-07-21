#ifndef MX_SHAPE
#define MX_SHAPE

// MxShape: abstract base class
// a shape is a geometric object defined by a function which returns a 
// positive value for points inside the shape, zero for points on the surface,
// and negative values for points outside the shape. 
//
//#include "MxConfig.h"
#include "Epetra_ConfigDefs.h"

#include <utility>
#include <string>

#include "MxUtil.hpp"
#include "MxDimVector.hpp"
#include "MxDimMatrix.hpp"
#include "MxGrid.h"

#include "Teuchos_XMLObject.hpp"


template<typename T, size_t DIM> class MxDimVector;

template<size_t DIM>
class MxShape {
  public:
    typedef std::pair<MxDimVector<double, DIM>, MxDimVector<double, DIM> > DimVecDblPair;

    virtual ~MxShape() {}

    virtual MxShape<DIM> * clone() const = 0;

    // 3D vecs can be initialized from 1D/2D vecs. 
    // For non-origin pivot point, the rotation is performed by a 
    // translation-rotation-translation operation.
    virtual void rotate(const MxVecD3 & axis, double angle, const MxVecD3 & pivot);

    // this is always a rotation of the shape about the given axis through the point 
    // to which the shape has been translated by means of this->translate(...)
    virtual void rotate(const MxVecD3 & axis, double angle);

    virtual void rotate(const Teuchos::XMLObject & node);

    // 3D vecs can be initialized from 1D/2D vecs 
    virtual void scale(const MxVecD3 & magnitudes);

    virtual void scale(const MxVecD3 & magnitudes, const MxVecD3 & origin);

    virtual void scale(const Teuchos::XMLObject & node);

    virtual void translate(const MxVecD3 & v);

    virtual void translate(const Teuchos::XMLObject & node);

    virtual void reflect(const MxVecD3 & normal, const MxVecD3 & pointInPlane = 0);

    virtual void reflect(const Teuchos::XMLObject & node);

    virtual void invert() {sign *= -1.0;}

    virtual void transforms(const Teuchos::XMLObject & node);

    virtual std::string getName() const {return name;}

    // 
    void applyTransforms(MxShape<DIM> & shape) const;

    // This is the defining function for the shape. It applies the Galilean transformations
    // that are general to
    // all shapes, then calls the purely virtual protected function '_func', giving as its
    // argument the transformed input coordinate. '_func' is
    // implemented by concrete shape subclasses.
    virtual double func(const MxDimVector<double, DIM> & p) const;

    // also return the 
    virtual double func(const MxDimVector<double, DIM> & p, int & indx) const;

    virtual MxDimVector<double, DIM> gradFunc(const MxDimVector<double, DIM> & p) const;

    // simply calls gradFunc, then returns the normalized result
    virtual MxDimVector<double, DIM> normal(MxDimVector<double, DIM> p) const;

    virtual int getNumSubShapes() const {return 1;}

    virtual MxShape<DIM> * getSubShape(int indx) const {return this->clone();}

    virtual bool hasGradFunc() const = 0;

    virtual MxDimVector<double, DIM> boundaryPoint(const MxDimVector<double, DIM> & p1, const MxDimVector<double, DIM> & p2) const;

    virtual DimVecDblPair boundingBox() const = 0;

    virtual void save(const MxGrid<DIM> & aGrid, const char * dSetName = "data", const char * fileName = 0) const;

  protected:

    // This function describes the untransformed shape. To be implemented by concrete
    // MxShape subclasses. Subclassed shapes should always be 3D, since they can be
    // transformed arbitrarily by MxShape transform functions.
    virtual double _func(const MxVecD3 & p) const = 0;

    virtual double _func(const MxVecD3 & p, int & indx) const {indx = 0; return this->_func(p);}

    virtual MxVecD3 _gradFunc(const MxVecD3 & p) const = 0;
    
    std::string name;

    double sign;

    MxDimVector<double, 3> b, Ainvb;

    MxDimMatrix<double, 3> A, Ainv;

    void init(const std::string & theName);

    void init(const Teuchos::XMLObject & node);
};

template<size_t DIM>
inline
void MxShape<DIM>::init(const std::string & theName) {
  name = theName;
  sign = 1.0;
  A = MxUtil::eye<double, 3>();
  Ainv = MxUtil::eye<double, 3>();
  b = 0;
  Ainvb = 0;
}

template<size_t DIM>
inline
void MxShape<DIM>::init(const Teuchos::XMLObject & node) {
  name = MxUtil::XML::getAttr("name", node);
  sign = 1.0;
  A = MxUtil::eye<double, 3>();
  Ainv = MxUtil::eye<double, 3>();
  b = 0;
  Ainvb = 0;
  transforms(node);
}

template<size_t DIM>
inline
double MxShape<DIM>::func(const MxDimVector<double, DIM> & p) const {
  MxVecD3 pp(p);
  return sign * this->_func(Ainv * pp - Ainvb);
}

template<size_t DIM>
inline
double MxShape<DIM>::func(const MxDimVector<double, DIM> & p, int & indx) const {
  MxVecD3 pp(p);
  //double f = sign * this->_func(Ainv * pp - Ainvb, indx);
  //std::cout << name << ": f = " << f << "\n";
  //return f;
  return sign * this->_func(Ainv * pp - Ainvb, indx);
}

template<size_t DIM>
inline
MxDimVector<double, DIM> MxShape<DIM>::gradFunc(const MxDimVector<double, DIM> & p) const {
  MxDimVector<double, 3> pp(p);
  return sign * MxDimVector<double, DIM>((this->_gradFunc(Ainv * pp - Ainvb)) * Ainv);
}

template<size_t DIM>
inline
MxDimVector<double, DIM> MxShape<DIM>::normal(MxDimVector<double, DIM> p) const {
  MxDimVector<double, DIM> res(this->gradFunc(p));
  return res / res.norm();
}


#endif
