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

# Sources of Core/Geometry classes

SET(Core_Geometry_SRCS
  BBox.cc
  CompGeom.cc
  Grid.cc
  Point.cc
  Plane.cc
  Polygon.cc
  Quaternion.cc 
  Ray.cc
  Transform.cc
  Vector.cc
)

SET(Core_Geometry_LIBS 
  Core_Persistent
  Core_Exceptions
  Core_Util
  Core_Math 
  ${M_LIBRARY}
  ${SLIVR_LIBRARY}
)

IF(SCI_TEEM_LIBRARY)
  # only use teem (and Tensor) if it's available
  #   This is mainly for non-SCIRun apps who don't use teem.
  #   Since SCIRun requires teem, this won't break it.  - BJW
  SET(Core_Geometry_SRCS
    ${Core_Geometry_SRCS}
    Tensor.cc
  )
  SET( Core_Geometry_LIBS ${Core_Geometry_LIBS} ${SCI_TEEM_LIBRARY} ${SCI_PNG_LIBRARY} ${SCI_ZLIB_LIBRARY})
ENDIF(SCI_TEEM_LIBRARY)

SCIRUN_ADD_LIBRARY(Core_Geometry ${Core_Geometry_SRCS})
TARGET_LINK_LIBRARIES(Core_Geometry ${Core_Geometry_LIBS})

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Core_Geometry)
ENDIF(BUILD_SHARED_LIBS)
