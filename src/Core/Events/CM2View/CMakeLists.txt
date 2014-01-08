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

# Sources of Core/Events/CM2View classes

SET(Core_Events_CM2View_SRCS
  CM2CursorTool.cc
  CM2ScaleTool.cc
  CM2TranslateTool.cc
  CM2UndoTool.cc
  CM2View.cc
  CM2WidgetTransformTool.cc
)

SCIRUN_ADD_LIBRARY(Core_Events_CM2View ${Core_Events_CM2View_SRCS})

TARGET_LINK_LIBRARIES( Core_Events_CM2View
  Core_Events
  Core_Volume
  ${GL_LIBRARY}
)

IF(WIN32)
  TARGET_LINK_LIBRARIES( Core_Events_CM2View 
                         ${GL_LIBRARY} )
ENDIF(WIN32)


IF(APPLE_LEOPARD)
  SET_TARGET_PROPERTIES(Core_Events_CM2View PROPERTIES LINK_FLAGS  "${GL_LINK}")
ENDIF(APPLE_LEOPARD)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Core_Events_CM2View)
ENDIF(BUILD_SHARED_LIBS)

