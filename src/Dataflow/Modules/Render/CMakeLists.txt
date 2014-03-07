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

# CMakeLists.txt for Dataflow/Modules/Render

SET(Dataflow_Modules_Render_SRCS
  Ball.cc
  BallAux.cc
  BallMath.cc
  EditPath.cc
  OpenGL.cc
  SynchronizeGeometry.cc
  ViewWindow.cc
  ViewSlices.cc
  ViewScene.cc
  glMath.cc
)

IF(BUILD_BIOMESH3D_REMOTE_SUPPORT)
	INCLUDE_DIRECTORIES(${MEDIA_COMMUNICATION_INCLUDE})
	INCLUDE_DIRECTORIES(${SCI_MEDIAMESSAGE_INCLUDE})
	INCLUDE_DIRECTORIES(${MEDIA_MESSAGE_INCLUDE})
	INCLUDE_DIRECTORIES(${BIOMESH3D_MESSAGE_INCLUDE})
	INCLUDE_DIRECTORIES(${JRTP_INCLUDE})
	INCLUDE_DIRECTORIES(${JTHREADS_INCLUDE})
	INCLUDE_DIRECTORIES(${FFMPEG_INCLUDE})
	INCLUDE_DIRECTORIES(${SCI_CRYPTO_INCLUDE})
	INCLUDE_DIRECTORIES(${SCI_CRYPTO_INCLUDE}/../include)
ENDIF(BUILD_BIOMESH3D_REMOTE_SUPPORT)

FOREACH(definition ${SCI_TK_DEFINITIONS})
  ADD_DEFINITIONS("-D${definition}")
ENDFOREACH(definition ${SCI_TK_DEFINITIONS})

SCIRUN_ADD_LIBRARY(Dataflow_Modules_Render ${Dataflow_Modules_Render_SRCS})

IF ( BUILD_BIOMESH3D_REMOTE_SUPPORT )
  SET( BIOMESH_REMOTE_LIBS
    ${SCI_BOOST_LIBRARY}
    ${JPEG_LIBRARY}
    ${SCI_MEDIAMESSAGE_LIBRARY}
    ${MEDIA_COMMUNICATION_LIBRARY}
    ${MEDIA_MESSAGE_LIBRARY}
    ${BIOMESH3D_MESSAGE_LIBRARY}
    ${JRTP_LIBRARY}
    ${JTHREADS_LIBRARY}
)
ENDIF( BUILD_BIOMESH3D_REMOTE_SUPPORT )

TARGET_LINK_LIBRARIES(Dataflow_Modules_Render
  Core_Algorithms_Fields
  Core_Basis
  Core_Datatypes
  Core_Events
  Core_Exceptions
  Core_Geometry
  Core_Geom
  Core_GeomInterface
  Dataflow_GuiInterface
  Core_Persistent
  Core_Thread
  Dataflow_TkExtensions
  Core_Util
  Core_Volume
  Dataflow_Comm
  Dataflow_Network
  Dataflow_Widgets
  ${TK_LIBRARY}
  ${GL_LIBRARY}
  ${M_LIBRARY}
  ${SCI_TEEM_LIBRARY}
  ${X_LIBRARY}
  ${BIOMESH_REMOTE_LIBS}
)

IF (BUILD_BIOMESH3D_REMOTE_SUPPORT)
	INCLUDE_DIRECTORIES(${JRTP_INCLUDE_DIR})
	INCLUDE_DIRECTORIES(${JTHREAD_INCLUDE_DIR})
ENDIF(BUILD_BIOMESH3D_REMOTE_SUPPORT)

IF(APPLE_LEOPARD)
  SET_TARGET_PROPERTIES(Dataflow_Modules_Render PROPERTIES LINK_FLAGS  "${GL_LINK} ${X_LINK}")
ENDIF(APPLE_LEOPARD)
IF(UNIX AND NOT WIN32)
    IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
      SET_TARGET_PROPERTIES(Dataflow_Modules_Render PROPERTIES COMPILE_FLAGS "-fPIC")
    ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
ENDIF(UNIX AND NOT WIN32)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Dataflow_Modules_Render)
ENDIF(BUILD_SHARED_LIBS)
