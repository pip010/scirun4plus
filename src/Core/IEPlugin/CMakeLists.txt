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

# Sources of Core/ImportExport classes

# NOTE: source files are all in sub directories ColorMap, Matrix, Field,
# and Nrrd
# 

SET(Core_IEPlugin_SRCS
  BYUSurfField_Plugin.cc
  CVRTIFileToMatrix_Plugin.cc
  CurveField_Plugin.cc
  DIF_Plugin.cc
  EcgsimFileToMatrix_Plugin.cc
  EcgsimFileToTriSurf_Plugin.cc
  HexVolField_Plugin.cc
  MRC_Plugin.cc
  MTextFileToTriSurf_Plugin.cc
  MatlabFiles_Plugin.cc
  NrrdField_Plugin.cc
  NrrdToMatrix_Plugin.cc
  ObjToField_Plugin.cc
  PointCloudField_Plugin.cc
  PVWave_Plugin.cc
  QuadSurfField_Plugin.cc
  SimpleTextFileToMatrix_Plugin.cc
  StructCurveField_Plugin.cc
  StructHexVolField_Plugin.cc
  StructQuadSurfField_Plugin.cc
  TriSurfField_Plugin.cc
  TetVolField_Plugin.cc
)


SCIRUN_ADD_LIBRARY(Core_IEPlugin ${Core_IEPlugin_SRCS})

TARGET_LINK_LIBRARIES(Core_IEPlugin
  Core_Persistent 
  Core_Exceptions 
  Core_Thread 
  Core_Geometry 
  Core_Math 
  Core_Util
  Core_Matlab 
  Core_Datatypes
  Core_ImportExport
  Core_Algorithms_DataIO
  Core_Algorithms_Converter 
  Core_Algorithms_Util
)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Core_IEPlugin)
ENDIF(BUILD_SHARED_LIBS)
