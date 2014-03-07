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

PROJECT(SLIVR)

IF(UNIX AND NOT APPLE AND NOT CYGWIN)
  SET(ENABLE_PIC ON)
ENDIF(UNIX AND NOT APPLE AND NOT CYGWIN)

IF (UNIX OR BUILD_STATIC_THIRDPARTY)
  ADD_DEFINITIONS(-DGLEW_STATIC)
ENDIF (UNIX OR BUILD_STATIC_THIRDPARTY)

# -- Source files the lib depends upon. --
# none of the following depend on GLEW
SET(SLIVR_SRCS
  BBox.cc
  Plane.cc
  Point.cc
  Ray.cc
  Transform.cc
  Vector.cc
  Color.cc
  ColorMap.cc
  ColorMap2.cc
  CM2Shader.cc
  CM2Widget.cc
  ShaderProgramARB.cc
  SliceRenderer.cc
  SurfaceRenderer.cc
  Texture.cc
  TextureBrick.cc
  TextureRenderer.cc
  VideoCardInfo.c
  VolShader.cc
  VolumeRenderer.cc
  )


SET(SCI_SLIVR_LIBRARY slivr CACHE INTERNAL "Slivr library name." FORCE)
SET(SCI_SLIVR_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/.. CACHE INTERNAL "Slivr include directories." FORCE)

INCLUDE_DIRECTORIES(
                    ${SCI_SLIVR_INCLUDE}
                    ${SCI_GLEW_INCLUDE}
                    ${SCI_TEEM_INCLUDE}
                   )

IF(UNIX OR BUILD_STATIC_THIRDPARTY)
  ADD_LIBRARY(${SCI_SLIVR_LIBRARY} STATIC ${SLIVR_SRCS})
ELSE(UNIX OR BUILD_STATIC_THIRDPARTY)
  ADD_LIBRARY(${SCI_SLIVR_LIBRARY} ${SLIVR_SRCS})
ENDIF(UNIX OR BUILD_STATIC_THIRDPARTY)

TARGET_LINK_LIBRARIES(${SCI_SLIVR_LIBRARY}
                      ${SCI_GLEW_LIBRARY}
                      ${SCI_TEEM_LIBRARY}
                      ${GL_LIBRARY}
                    )

IF(ENABLE_PIC)
  SET_TARGET_PROPERTIES(${SCI_SLIVR_LIBRARY}
                        PROPERTIES
                        COMPILE_FLAGS "-fPIC")
ENDIF(ENABLE_PIC)
