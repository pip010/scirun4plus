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

# Sources of Core/Events classes

SET(Core_Events_SRCS
  BaseEvent.cc
  BaseTool.cc
  CursorChangeEvent.cc
  DataManager.cc
  EventManager.cc
  OpenGLBase.cc
  SceneGraphEvent.cc
  ToolManager.cc	
  OSXEventSpawner.cc
)

IF(WIN32)
  SET(Core_Events_SRCS ${Core_Events_SRCS} Win32EventSpawner.cc)
ELSE(WIN32)
  SET(Core_Events_SRCS ${Core_Events_SRCS} X11EventSpawner.cc)
ENDIF(WIN32)

# DV Fix should do subdirectory tools instead of adding tools files above

ADD_SUBDIRECTORY(CM2View)
ADD_SUBDIRECTORY(OGLView)

SCIRUN_ADD_LIBRARY(Core_Events ${Core_Events_SRCS})

TARGET_LINK_LIBRARIES( Core_Events
  Core_Exceptions
  Core_Thread
  Core_Math
  Core_Geom
  Core_Util
  Core_Geometry
  Core_Datatypes
  ${THREAD_LIBRARY}
  ${GL_LIBRARY}
  ${SCI_PNG_LIBRARY}
  ${SCI_ZLIB_LIBRARY}
  ${SCI_TEEM_LIBRARY}
  ${SCI_FREETYPE_LIBRARY}
  ${M_LIBRARY}
  ${SCI_GLEW_LIBRARY}
)


IF(APPLE)
  IF(NOT APPLE_LEOPARD OR NOT WITH_X11)
    TARGET_LINK_LIBRARIES( Core_Events "-framework Carbon")    
  ELSE(NOT APPLE_LEOPARD OR NOT WITH_X11)
    SET_TARGET_PROPERTIES(Core_Events PROPERTIES LINK_FLAGS  "${GL_LINK} ${X_LINK}")
  ENDIF(NOT APPLE_LEOPARD OR NOT WITH_X11)
ENDIF(APPLE)

IF(WIN32)
  TARGET_LINK_LIBRARIES( Core_Events )
ENDIF(WIN32)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Core_Events)
ENDIF(BUILD_SHARED_LIBS)

