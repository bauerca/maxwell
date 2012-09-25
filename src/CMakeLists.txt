
# utils
set(headers
  MxTypes.h
  MxUtil.hpp
  MxDimMatrix.hpp
  MxDimVector.hpp
  MxDynMatrix.hpp
)

set(sources
  MxUtil.cpp
)

# fundamental shape geometries
set(headers ${headers}
  MxShape.hpp
  MxCylinder.hpp
  MxSphere.hpp
  MxTorus.hpp
  MxCone.hpp
  MxSlab.hpp
  MxEllipsoid.hpp
  MxHalfSpace.hpp
  MxShapeIntersection.hpp 
  MxShapeMirror.hpp      
  MxShapeRepeat.hpp      
  MxShapeSubtract.hpp    
  MxShapeUnion.hpp       
  MxShapeDifference.hpp
  MxShapeDefs.h
)

set(sources ${sources}
  MxShape.cpp
)

# cartesian grid geometric objects
set(headers ${headers}
  MxPoint.hpp
  MxPolytope.hpp
  MxConvexPolyhedron.h
  MxConvexPolygon.h
  MxCartRect.hpp
  MxCartSeg.hpp
  MxCartBox.hpp
  MxSegment.h
)

set(sources ${sources}
  MxPolytope.cpp
  MxCartBox.cpp
  MxSegment.cpp
  MxConvexPolygon.cpp
  MxConvexPolyhedron.cpp
)

# grid stuff
set(headers ${headers}
  MxGrid.h
  MxGridDomain.h
  MxGridDomainIter.hpp
)

set(sources ${sources}
  MxGrid.cpp
  MxGridDomain.cpp
)

# fields
set(headers ${headers}
  MxGridField.hpp
  MxGridFieldIter.hpp
  MxYeeElecFieldBase.h
  MxYeeFitEField.h
  MxYeeFitDField.h
  MxYeeMagFieldBase.h
  MxYeeFitBField.h
  MxYeeFitHField.h
  MxYeePsiField.h
)

set(sources ${sources}
  MxGridField.cpp     
  MxYeeElecFieldBase.cpp
  MxYeeFitEField.cpp
  MxYeeFitDField.cpp
  MxYeeMagFieldBase.cpp
  MxYeeFitBField.cpp
  MxYeeFitHField.cpp
  MxYeePsiField.cpp
)

# linear algebra
set(headers ${headers}
  MxMap.hpp
  MxMultiVector.hpp
  MxAnasaziMV.hpp
  MxCrsMatrix.hpp
)

set(sources ${sources}
  MxMap.cpp
  MxMultiVector.cpp
  MxAnasaziMV.cpp
  MxCrsMatrix.cpp
)

# operators
set(headers ${headers}
  MxOperator.hpp
  MxGridFieldInterpolator.h
  MxYeeDeyMittraCurlB.h
  MxYeeDeyMittraCurlE.h
  MxYeeDeyMittraFracs.h
  MxYeeDeyMittraGradPsi.h
  MxYeeDeyMittraDivB.h
  MxYeeFitInvEps.h
  MxYeeFitMu.h
  #MxYeeMitInvEps.h
  MxGridFieldFuncOp.h
  MxEikxOp.h
  MxEMOps.h
  MxMagWaveOp.h
)

set(sources ${sources}
  MxGridFieldInterpolator.cpp
  MxYeeDeyMittraCurlB.cpp
  MxYeeDeyMittraCurlE.cpp
  MxYeeDeyMittraFracs.cpp
  MxYeeDeyMittraGradPsi.cpp
  MxYeeDeyMittraDivB.cpp
  MxYeeFitInvEps.cpp
  MxYeeFitMu.cpp
  #MxYeeMitInvEps.cpp
  MxGridFieldFuncOp.cpp
  MxEikxOp.cpp
  MxEMOps.cpp
  MxMagWaveOp.cpp
)

# high-level objects
set(headers ${headers}
  MxDielectric.hpp
  MxMu.hpp
  MxPML.h
  MxPointCloud.h
  MxEMSim.h
  #MxEMSimHierarchy.h
  #MxGeoMultigridPrec.h
  MxIO.h
  MxSolver.h
  MxProblem.h
)

set(sources ${sources}
  MxDielectric.cpp
  MxMu.cpp
  MxPML.cpp
  MxPointCloud.cpp
  MxEMSim.cpp
  #MxEMSimHierarchy.cpp
  #MxGeoMultigridPrec.cpp
  MxIO.cpp
  MxSolver.cpp
  MxProblem.cpp
)


add_library(mxwll ${headers} ${sources})

add_executable(maxwell maxwell.cpp)
#add_executable(maxwell maxwell.cpp ${headers} ${sources})

add_dependencies(maxwell MxHDF5 MxTrilinos)

#target_link_libraries(maxwell ${Trilinos_LIBRARIES} hdf5 blas lapack z m)
target_link_libraries(maxwell mxwll
  ${MAXWELL_EXTPROJ_LINK_LIBS} blas lapack z m)


install(TARGETS maxwell RUNTIME DESTINATION bin)
