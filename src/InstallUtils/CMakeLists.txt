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

IF (WITH_TETGEN)
  SET(SCIRUN_BUILD_DIST "SCIRun_${SCIRUN_VERSION_MAJOR}.${SCIRUN_VERSION_MINOR}")
ELSE(WITH_TETGEN)  
  SET(SCIRUN_BUILD_DIST "SCIRunNoMeshing_${SCIRUN_VERSION_MAJOR}.${SCIRUN_VERSION_MINOR}")
ENDIF(WITH_TETGEN)

IF(WIN32)
  ADD_SUBDIRECTORY(windows)
ENDIF(WIN32)

IF(APPLE)
  FIND_PROGRAM(PLATYPUS_EXECUTABLE platypus)

  SET(SCIRUN_ICON_FILE "${CMAKE_SOURCE_DIR}/Applications/SCIRun/SCIRun-icon.icns")
  SET(SCIRUN_STARTUP_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/startup.sh")
  SET(SCIRUN_APP_NAME SCIRun.app)

  SET(SEG3D_APP_NAME Seg3D.app)

  SET(BIOMESH3D_APP_NAME BioMesh3D.app)
  SET(BIOMESH3D_STARTUP_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/BioMesh3D.sh")

  #SET(BIOTENSOR_APP_NAME BioTensor.app)
  #SET(BIOMESH3D_STARTUP_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/BioTensor.sh")

  IF (CMAKE_SIZEOF_VOID_P MATCHES 8)
    SET(SCIRUN_BUILD_TYPE "osx64")
  ELSE (CMAKE_SIZEOF_VOID_P MATCHES 8)
    SET(SCIRUN_BUILD_TYPE "osx32")
  ENDIF (CMAKE_SIZEOF_VOID_P MATCHES 8)

  IF(ITK_FOUND AND ITK_BUILD_SHARED_LIBS MATCHES "SHARED")
    SET(DISABLE_ITK "")
  ELSE(ITK_FOUND AND ITK_BUILD_SHARED_LIBS MATCHES "SHARED")
    SET(DISABLE_ITK "#")
  ENDIF(ITK_FOUND AND ITK_BUILD_SHARED_LIBS MATCHES "SHARED")

  SET(SEG3D_BUILD_DIST "Seg3D_${SEG3D_VERSION_MAJOR}.${SEG3D_VERSION_MINOR}")

  CONFIGURE_FILE(buildApp.sh.in ${CMAKE_BINARY_DIR}/buildApp.sh @ONLY)

  CONFIGURE_FILE(buildSeg3D.sh.in ${CMAKE_BINARY_DIR}/buildSeg3D.sh @ONLY)

  CONFIGURE_FILE(BioTensor.sh ${CMAKE_BINARY_DIR}/InstallUtils/BioTensor @ONLY)
ENDIF(APPLE)
