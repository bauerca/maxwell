#ifndef MX_GRID_FIELD
#define MX_GRID_FIELD

#include <vector>
#include <string>
#include <map>

#include "Epetra_ConfigDefs.h"

#include "MxTypes.h"
#include "MxDimVector.hpp"
#include "MxPolytope.hpp"
#include "MxShape.hpp"
#include "MxGrid.h"
#include "MxGridDomainIter.hpp"
#include "MxPointCloud.h"

#include "MxVector.hpp"

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

template<size_t> class MxGridFieldIter;

/**
A grid field is defined on a Cartesian grid. It has an arbitrary
number of field components (e.g. an electric field would have 3).
Each field component is given a Polytope representation (could just
be a point for the finite difference method, but could also be a
line or polygon if using the finite integration technique) and a
cell coordinate (cell-relative location of the centroid of the
field component's associated polytope).

The grid field can hold representations of shape objects as given by
the intersection of the shape object with each component of the grid
field. Specifically, the grid field holds a fraction (0-1) for each
component that gives the percentage of the "volume" of the associated
polytope \em inside the shape.

The grid field is responsible for generating a Map which holds all of
the \em relevant component indices (an irrelevant component would
always have a value equal to zero; e.g. a
component for an electric field that lies within a perfect conductor).
The Map is used in the construction of matrices and vectors that are
passed on to linear algebra libraries like Trilinos.
*/
template<size_t DIM>
class MxGridField {
  public:
    //MxGridField(const MxGrid<DIM> * aGrid);

    //virtual size_t addComp(const MxPolytope<DIM> & aPtope, const MxDimVector<double, DIM> & aCellCoord);

    virtual ~MxGridField() {}

    virtual std::string getName() const {return fieldName;}

    virtual const MxGrid<DIM> & getGrid() const {return *grid;}

    virtual const MxPolytope<DIM> & getCompPolytope(size_t comp) const = 0;

    /**
    Grid-cell-relative position of field component. This position should
    \em never coincide with an upper bound of a grid cell.
    */
    virtual const MxDimVector<double, DIM> & getCompCellCoord(size_t comp) const = 0;

    virtual void addShapeRep(const MxShape<DIM> & aShape, std::string shapeName, int numGuard = 1, bool setRegion = false);

    virtual void addShapeRep(const std::vector<double> & compFracs, std::string shapeName, int numGuard, bool setRegion = false) {;}

    virtual bool hasShapeRep(std::string shapeName) const;

    //virtual void setRegion(const MxShape<DIM> & aShape);
    //virtual void setRegion(size_t shapeIndx);
    virtual void setRegion(std::string shapeName);

    virtual void setCompBCs(size_t comp,
      MxDimVector<MxBCType, DIM> lowerBCs,
      MxDimVector<MxBCType, DIM> upperBCs);

    /**
     * Higher-level set boundary conditions function for setting bcs on all components
     * simultaneously (e.g. PEC/PMC for electromagnetic fields)
     */
    //virtual void setBCs(MxDimVector<MxBCType, DIM> lowerBCs,
    //    MxDimVector<MxBCType, DIM> upperBCs) = 0;

    virtual void setPhaseShifts(MxDimVector<double, DIM> phaseShifts);

    virtual bool useCompInMap(size_t comp, MxDimVector<int, DIM> cell) const;

    virtual bool compInMap(size_t comp, MxDimVector<int, DIM> cell) const;

    virtual size_t getNumComps() const {return numComps;}

    /**
    Build the Map.
    */
    virtual void setMap();

    /**
    Set the Map for this field to the Map of another field. This is useful
    when fields are basically the same thing, such as the E,D and B,H
    fields in Maxwell's equations.
    */
    virtual void setMap(const MxGridField<DIM> & partnerField);

    virtual RCP<MxMap> getMap() const {return mMap;}

    /**
    Get the unique index for the given component. If a component is
    provided that is outside the bounds of the simulation, the
    component is first "interiorized" (using
    MxGridField<size_t>::getInteriorComp()) according to the boundary
    conditions set on the component (e.g. Dirichlet/Neumann BCs will
    \em reflect the component back to the interior while periodic
    BCs \em wrap the component). The result of
    MxGridField<size_t>::getCompFactor() using the same component and
    cell indices is the coefficient that should be applied to the
    interior component value so that the boundary conditions are
    satisfied.
    
    For example, Dirichlet BCs that set a component to zero on the
    boundary would cause an out-of-bounds component to be reflected
    back to the interior with a coefficient of -1.

    Beware, reflected components may have different cell coordinates
    if the original cell coordinates are not centered or located on
    a cell boundary in any direction. These cases are not treated since
    the foreseeable future will use the Yee layout.
    */
    virtual size_t globCompIndx(size_t comp, const MxDimVector<int, DIM> & cell) const {
      return comp + numComps * grid->cellToGlobalIndx(
          this->getInteriorComp(comp, cell));
    }

    /**
    Get the coefficient for the component index returned by
    MxGridField<size_t>::globCompIndx(). The coefficient will be
    equal to 1 unless the comp/cell combination refers to a 
    component that is outside the simulation domain. See documentation
    for MxGridField<size_t>::globCompIndx() for many more details.
    */
    virtual MxComplex getCompFactor(size_t comp, MxDimVector<int, DIM> cell) const;

    /**
    Return the cell containing the "interiorized" component. That is,
    given a comp/cell combination that refers to a component outside
    the simulation bounds, this function returns the appropriate
    in-bounds containing cell; a component is interiorized according
    to simulation boundary conditions.

    Only the cell is returned since the field component remains the
    same. For nonstandard field component cell coordinates, this can
    raise problems. See MxGridField<size_t>::globCompIndx() for more
    info.
    */
    virtual MxDimVector<int, DIM> getInteriorComp(size_t comp,
      MxDimVector<int, DIM> cell) const;

    virtual MxDimVector<double, DIM> getCompCoord(size_t comp, const MxDimVector<int, DIM> & cell) const {
      return grid->nodeCoord(cell) + compCellCoords[comp];
    }

    /**

    */
    virtual double calcCompFrac(size_t comp, MxDimVector<int, DIM> cell,
        const MxShape<DIM> & aShape) const {
      MxDimVector<int, DIM> newCell(this->getInteriorComp(comp, cell));
      return compPtopes[comp]->volumeFraction(aShape,
          grid->nodeCoord(newCell) + compCellCoords[comp]);
    }

    //virtual double getCompFrac(size_t comp, const MxDimVector<int, DIM> & cell, size_t shapeIndx) const;
    virtual double getCompFrac(size_t comp, const MxDimVector<int, DIM> & cell, std::string shapeName) const;

    //virtual Epetra_Vector getAllCompFracs(int shapeIndx) const;
    virtual RCP<MxVector<double> > getAllCompFracs(std::string shapeName) const;

    template<typename Scalar>
    static RCP<MxMultiVector<Scalar> > uniformFields(
      MxGridField<DIM> const & field);

    virtual RCP<MxMultiVector<double> > coords() const;

    friend class MxGridFieldIter<DIM>;

  protected:
    std::string fieldName;

    const MxGrid<DIM> * grid;

    size_t numComps;

    //Teuchos::RCP<MxMap> map;
    RCP<MxMap> mMap;

    bool regionSet;

    //int region;
    std::string regionName;

    bool mapSet;

    std::vector<Teuchos::RCP<MxPolytope<DIM> > > compPtopes;

    std::vector<MxDimVector<double, DIM> > compCellCoords;

    std::vector<MxDimVector<MxBCType, DIM> > mCompLowerBCs;

    std::vector<MxDimVector<MxBCType, DIM> > mCompUpperBCs;

    MxDimVector<MxComplex, DIM> mPhaseFactors;

    //std::vector<std::vector<double> > shapeReps;
    std::map<std::string, std::vector<double> > shapeReps;
    //std::map<std::string, Teuchos::RCP<Epetra_MultiVector> > shapeReps;
    //std::map<std::string, Teuchos::RCP<Epetra_Vector> > shapeReps;
    //std::map<std::string, MxGridFieldData<DIM> > shapeReps;

    //std::vector<std::string> shapeRepNames;

    //std::vector<MxGridDomain<DIM> > shapeRepsDomains;
    std::map<std::string, MxGridDomain<DIM> > shapeRepsDomains;

    Teuchos::RCP<Epetra_MultiVector> fields;

    void initBCs();

};

//template<size_t DIM>
//MxGridField<DIM>::MxGridField(const MxGrid<DIM> * aGrid) : grid(aGrid), domain(aGrid->getDomain()), mapSet(false), regionSet(false) {}


template<size_t DIM>
inline
double MxGridField<DIM>::getCompFrac(size_t comp, const MxDimVector<int, DIM> & cell, std::string shapeName) const {
  //if (shapeReps.size() == 0 or regionSet == false) return 1;
  //else if (shapeIndx < 0 or shapeIndx > shapeReps.size() - 1) {
  //  std::cout << "MxGridField<" << DIM << ">::getCompVolFrac: shapeIndx (" 
  //            << shapeIndx << ") out of range (0 - " << shapeReps.size() - 1 << ").\n";
  //  throw 1;
  //}
  //else {
  //  MxGridDomain<DIM> domain(grid->getGridDomain(shapeRepsGuardCells[shapeIndx]));
  //  return shapeReps[shapeIndx][comp + numComps * domain.cellToFullIndx(cell)];
  //}
  
  if (hasShapeRep(shapeName))
    return shapeReps.find(shapeName)->second[comp + numComps * shapeRepsDomains.find(shapeName)->second.cellToFullIndx(cell)];
  else
    return 1.0;
}

template<size_t DIM>
inline
bool MxGridField<DIM>::hasShapeRep(std::string shapeName) const {
  return shapeReps.find(shapeName) != shapeReps.end();
}


#endif
