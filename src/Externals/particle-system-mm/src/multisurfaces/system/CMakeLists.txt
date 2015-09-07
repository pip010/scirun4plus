SET(system_SRCS
 mtxlib.cpp
 ParticleSystem.cxx
 SFSystem.cxx
 optimization/AdaptiveDtSurfaceDistribution.cxx
 optimization/ApproximateNeighborhoods.cxx
 optimization/GlobalSurfaceEnergy.cxx
 optimization/LODController.cxx
 optimization/Optimize.cxx
 optimization/SurfaceConstraint.cxx
 optimization/SurfacePopulationController.cxx
 particles/DynamicSurfacePoint.cxx
 particles/Point.cxx
 particles/SurfacePoint.cxx
 domain/DistanceField.cxx
 domain/DistanceFieldwSF.cxx
 domain/Domain.cxx
 domain/Energy.cxx
 domain/Intersection.cxx
 domain/IOScalarField.cxx
 domain/ScalarField.cxx
 domain/SimpleSurfaces.cxx
 domain/Surface.cxx
 domain/SurfaceParameters.cxx
)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}  ${PARTICLEUTIL_INCLUDES})

ADD_LIBRARY(system ${system_SRCS})
TARGET_LINK_LIBRARIES(system particleutil)
IF(WIN32)
  SET_TARGET_PROPERTIES(system
                        PROPERTIES
                          COMPILE_FLAGS "-DBUILD_system"
                          DEFINE_SYMBOL _USRDLL
                       )
ENDIF(WIN32)
