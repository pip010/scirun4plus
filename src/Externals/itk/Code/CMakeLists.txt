SUBDIRS(Common Numerics BasicFilters Algorithms SpatialObject)

#######################################
# SCI - IGNORE THESE ONES
#IF( ITK_USE_PATENTED )
#  SUBDIRS( Patented )
#ENDIF( ITK_USE_PATENTED )
#
IF( ITK_USE_REVIEW )
  SUBDIRS( Review )
ENDIF( ITK_USE_REVIEW )

## There should be no library dependances on IO, but IO may depend on other
## libraries, so build it last
SUBDIRS(IO)

SET(SCI_ITK_LIBRARY "ITKAlgorithms;ITKBasicFilters;ITKNumerics;ITKFEM;ITKSpatialObject" CACHE INTERNAL "ITK libraries library name." FORCE)
SET(SCI_ITK_LIBRARIES "ITKAlgorithms;ITKBasicFilters;ITKNumerics;ITKFEM;ITKSpatialObject" CACHE INTERNAL "ITK libraries library name." FORCE)
SET(SCI_ITK_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/Common;${CMAKE_CURRENT_SOURCE_DIR}/Numerics/FEM;${CMAKE_CURRENT_SOURCE_DIR}/Numerics;${CMAKE_CURRENT_SOURCE_DIR}/BasicFilters;${CMAKE_CURRENT_SOURCE_DIR}/Algorithms;${CMAKE_CURRENT_SOURCE_DIR}/SpatialObject;${CMAKE_CURRENT_SOURCE_DIR}/IO;${ITK_INCLUDE_DIRS_BUILD_TREE};${ITK_INCLUDE_DIRS_BUILD_TREE_CXX};${ITK_INCLUDE_DIRS_SYSTEM}" CACHE INTERNAL "ITK include directories" FORCE)