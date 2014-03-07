SET(optimization_SRCS
  AdaptiveDtSurfaceDistribution.cxx
  ApproximateNeighborhoods.cxx
  GlobalSurfaceEnergy.cxx
  LODController.cxx
  Optimize.cxx
  SurfaceConstraint.cxx
  SurfacePopulationController.cxx
  mtxlib.cpp
)

ADD_LIBRARY(optimization ${optimization_SRCS})

TARGET_LINK_LIBRARIES(optimization domain)

