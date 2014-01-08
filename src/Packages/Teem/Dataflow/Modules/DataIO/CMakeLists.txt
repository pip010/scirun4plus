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

# CMakeLists.txt for Packages/Teem/Dataflow/Modules/DataIO

SET(Packages_Teem_Dataflow_Modules_DataIO_SRCS
  ExportImage.cc
  ImportImage.cc
  ReadNrrd.cc
  WriteNrrd.cc
)

IF(ITK_FOUND)
  SET(Packages_Teem_Dataflow_Modules_DataIO_SRCS ${Packages_Teem_Dataflow_Modules_DataIO_SRCS}
  ImportNrrdFromAnalyze.cc
  ImportNrrdFromDicom.cc
)
ENDIF(ITK_FOUND)

SCIRUN_ADD_LIBRARY(Packages_Teem_Dataflow_Modules_DataIO ${Packages_Teem_Dataflow_Modules_DataIO_SRCS})

TARGET_LINK_LIBRARIES(Packages_Teem_Dataflow_Modules_DataIO
  Core_Datatypes
  Core_Algorithms_DataIO
  Core_Exceptions
  Core_Geom
  Core_Geometry
  Core_GeomInterface
  Dataflow_GuiInterface
  Core_ImportExport
  Core_Persistent
  Core_Thread
  Dataflow_TkExtensions
  Core_Util
  Dataflow_Network
  ${SCI_TEEM_LIBRARY}
  ${SCI_PNG_LIBRARY}
  ${SCI_ZLIB_LIBRARY}
)

IF(ITK_FOUND)
  TARGET_LINK_LIBRARIES(Packages_Teem_Dataflow_Modules_DataIO
  ITKIO
  ITKCommon
  ${GDCM_LIBRARY}
  )
ENDIF(ITK_FOUND)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Packages_Teem_Dataflow_Modules_DataIO)
ENDIF(BUILD_SHARED_LIBS)
