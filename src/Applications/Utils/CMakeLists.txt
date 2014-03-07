#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2009 Scientific Computing and Imaging Institute,
#  University of Utah.
# 
#  
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#

PROJECT(APPLICATIONS_UTILS)

SET(Utils_LIBS
  Core_Algorithms_DataIO
  Core_Algorithms_Fields
  Core_Algorithms_Util
  Core_Exceptions
  Core_Datatypes
  Core_Util
  Core_Persistent
  Core_Init
  Core_Thread
  Core_Geometry
  Core_Math
  Core_ImportExport
  Core_Basis
  Core_Geom
  ${SCI_TEEM_LIBRARY}
  ${SCI_FREETYPE_LIBRARY}
  ${SCI_LIBXML2_LIBRARY}
  ${SCI_PNG_LIBRARY}
  ${SCI_ZLIB_LIBRARY}
  ${M_LIBRARY}
  ${GL_LIBRARY}
  ${GLU_LIBRARY}
)

SET(Utils_SRCS 
  AddTri
  ApplyMappingMatrix
  BuildMappingMatrix
  ColumnMatrixToText
  ComputeTightenedLabels
  ConvertFieldToNrrd
  ConvertNrrdToField
  DenseMatrixToText
  ExtractIsosurface
  FieldBin1Test
  FieldTextToBin
  GenerateMPMData
  GenerateMedialAxisPoints
  JitterPointCloud
  JoinFields
  LabelTets
  MtoTriSurfField
  OrientFaces
  ObjToField
  PCA-example
  RemoveConnectedFaces
  RemoveFaces
  RemoveOrphanNodes
  SeedIntersections
  MapFieldDataFromSourceToDestination
  SparseRowMatrixToText
  ColorMapToText
  MA2PointCloudVectors
  VTKtoTriSurfField
  VTKtoHexVolField
  TextToColumnMatrix
  TextToDenseMatrix
  TextToColorMap
  TriSurfToOBJ
  TriSurfToTetgen
  TetgenToTetVol
  vff2nrrd
  nrrd2vff
  nrrd2vol
  PictToNrrd
  PtclToField
  FieldToPtcl
  VistaToNrrd
  gmvToPts
  geoToPts
  EGItoMat
  TextToCurve
  UnorientNrrdAndGetTransform
  TransformFieldWithTransform
)

IF(WITH_TETGEN)
  SET(Utils_SRCS ${Utils_SRCS} InterfaceWithTetGen)
  SET(Utils_LIBS ${Utils_LIBS} ${SCI_TETGEN_LIBRARY})
ENDIF(WITH_TETGEN)

FOREACH(util ${Utils_SRCS})
  
  ADD_EXECUTABLE(${util} ${SCIRUN_CORE_SOURCE_DIR}/Applications/Utils/${util}.cc)
  TARGET_LINK_LIBRARIES(${util} ${Utils_LIBS})
  
  #INSTALL(TARGETS ${util}
  #         RUNTIME DESTINATION Utils
  #       )
ENDFOREACH(util ${Utils_SRCS})



IF(APPLE_LEOPARD)
  SET_TARGET_PROPERTIES(ExtractIsosurface PROPERTIES LINK_FLAGS  "${GL_LINK} ${X_LINK}")
ENDIF(APPLE_LEOPARD)
