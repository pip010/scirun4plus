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

# CMakeLists.txt for Dataflow/Modules/Visualization

SET(Dataflow_Modules_Visualization_SRCS
  CreateScaleBar.cc
  CreateLightForViewer.cc
  ChooseColorMap.cc
  CreateAndEditColorMap.cc
  CreateAndEditColorMap2D.cc
  CreateViewerAxes.cc
  CreateViewerClockIcon.cc
  CreateStandardColorMaps.cc
  ExtractIsosurface.cc
  ConvertNrrdsToTexture.cc
  ConvertFieldsToTexture.cc
  RescaleColorMap.cc
  ShowColorMap.cc
  ShowFieldGlyphs.cc
  ShowField.cc
  ShowLinePath.cc
  ShowMatrix.cc
  ShowMeshBoundingBox.cc
  ShowPointPath.cc
  GenerateStreamLines.cc
  GenerateStreamLinesWithPlacementHeuristic.cc
  ShowTextureSlices.cc
  ShowTextureSurface.cc
  ShowTextureVolume.cc
  ShowString.cc
)

SCIRUN_ADD_LIBRARY(Dataflow_Modules_Visualization ${Dataflow_Modules_Visualization_SRCS})

TARGET_LINK_LIBRARIES(Dataflow_Modules_Visualization
  Core_Algorithms_Visualization
  Core_Algorithms_Fields
  Core_Algorithms_Math
  Core_Basis
  Core_Datatypes
  Core_Events
  Core_Events_CM2View
  Core_Exceptions
  Core_Geom
  Core_GeomInterface
  Core_Geometry
  Dataflow_GuiInterface
  Core_Persistent
  Core_Thread
  Dataflow_TkExtensions
  Core_Util
  Core_Volume
  Dataflow_Network
  Dataflow_Widgets
  ${SCI_TEEM_LIBRARY}
  ${FTGL_LIBRARY}
  ${TK_LIBRARY}
  ${GL_LIBRARY}
  ${M_LIBRARY}
)

IF(APPLE_LEOPARD)
  SET_TARGET_PROPERTIES(Dataflow_Modules_Visualization PROPERTIES LINK_FLAGS  "${GL_LINK}")
ENDIF(APPLE_LEOPARD)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Dataflow_Modules_Visualization)
ENDIF(BUILD_SHARED_LIBS)
