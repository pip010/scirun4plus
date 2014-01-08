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

# Sources of Core/Volume classes

SET(Core_Volume_SRCS
  CM2Widget.cc
  ColorMap2.cc
  SliceRenderer.cc
  SurfaceRenderer.cc
  Texture.cc
  TextureBrick.cc
  VolumeRenderer.cc
)


INCLUDE_DIRECTORIES(${SLIVR_INCLUDE_DIR})

SCIRUN_ADD_LIBRARY(Core_Volume ${Core_Volume_SRCS})


TARGET_LINK_LIBRARIES( Core_Volume 
  Core_Datatypes
  Core_Exceptions
  Core_Geom
  Core_Geometry
  Core_Persistent
  Core_Thread
  Core_Util
  ${GL_LIBRARY}
  ${M_LIBRARY}
  ${SCI_TEEM_LIBRARY}
  ${SLIVR_LIBRARY}
  ${SCI_GLEW_LIBRARY}
)

IF(APPLE_LEOPARD)
  SET_TARGET_PROPERTIES(Core_Volume PROPERTIES LINK_FLAGS  "${GL_LINK}")
ENDIF(APPLE_LEOPARD)

IF(WIN32 AND NOT CYGWIN)
  TARGET_LINK_LIBRARIES( Core_Volume Core_Math )
ENDIF(WIN32 AND NOT CYGWIN)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Core_Volume)
ENDIF(BUILD_SHARED_LIBS)
