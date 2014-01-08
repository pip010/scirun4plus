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

PROJECT(SLIVR_EXAMPLES)

# Force static builds
SET(BUILD_SHARED_LIBS OFF)

# Set which components of wx we need, as there can be different parts
# compiled.

SET(wxWidgets_INCLUDE_DIRS_NO_SYSTEM 1)
SET(wxWidgets_USE_LIBS core base gl)
FIND_PACKAGE(wxWidgets)

IF(NOT wxWidgets_FOUND)
  MESSAGE(FATAL_ERROR "Required wxWidget libraries not found")
ENDIF(NOT wxWidgets_FOUND)
INCLUDE(${wxWidgets_USE_FILE})


IF(APPLE)
  ADD_DEFINITIONS(-DEX_DATA_PATH="\\"./wxgladetest.app/data/\\"")
  INSTALL(FILES tooth.nhdr tooth.raw DESTINATION bin/examples/wxgladetest.app/data)
ELSE(APPLE)
  ADD_DEFINITIONS(-DEX_DATA_PATH="\\"../data/\\"")
  INSTALL(FILES tooth.nhdr tooth.raw DESTINATION bin/data)
ENDIF(APPLE)



# and for each of your dependant executable/library targets:
# MACOSX_BUNDLE : add this to make a MAC bundle.
ADD_EXECUTABLE(wxgladetest WIN32 MACOSX_BUNDLE
  main.cpp MainWindow.cpp OGLCanvas.cc)

TARGET_LINK_LIBRARIES(wxgladetest slivrV slivrG
  ${wxWidgets_LIBRARIES} ${TEEM_LIBRARIES} ${GLEW_LIBRARY})
IF(APPLE_LEOPARD)
  SET_TARGET_PROPERTIES(wxgladetest PROPERTIES LINK_FLAGS "${GL_LINK}")
ENDIF(APPLE_LEOPARD)

INCLUDE_DIRECTORIES(.)
INSTALL(TARGETS wxgladetest RUNTIME DESTINATION bin/examples)
