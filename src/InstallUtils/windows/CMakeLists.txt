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

# This configuration writes the setup executable to SCIRUN_CORE_BINARY_DIR.

FIND_PACKAGE(InnoSetup)

IF(INNO_SETUP_COMMAND)
  FILE(TO_NATIVE_PATH "${INNO_SETUP_COMMAND}" INNO_SETUP_COMMAND_NATIVE)
ELSE(INNO_SETUP_COMMAND)
  MESSAGE(STATUS "Warning: InnoSetup was not found. Creating a Windows executable installer will not work.")
ENDIF(INNO_SETUP_COMMAND)

SET(SHORT_COPYRIGHT "Copyright (c) 2009 Scientific Computing and Imaging Institute, University of Utah.")
SET(SCI_URL "http://software.sci.utah.edu")

FILE(TO_NATIVE_PATH "${SCIRUN_CORE_BINARY_DIR}" SCIRUN_BINARY_DIR_NATIVE)
STRING(REPLACE "/src" "" SCIRUN_TOP_DIR ${SCIRUN_CORE_SOURCE_DIR})
FILE(TO_NATIVE_PATH "${SCIRUN_CORE_SOURCE_DIR}" SCIRUN_SOURCE_DIR_NATIVE)
FILE(TO_NATIVE_PATH "${SCIRUN_TOP_DIR}" SCIRUN_TOP_DIR_NATIVE)
FILE(TO_NATIVE_PATH "${SCIRUN_THIRDPARTY_DIR}" SCIRUN_THIRDPARTY_DIR_NATIVE)
FILE(TO_NATIVE_PATH "${EXECUTABLE_OUTPUT_PATH}" EXECUTABLE_OUTPUT_PATH_NATIVE)

IF(CMAKE_GENERATOR MATCHES "Visual Studio 9")
  SET(VS_VERSION "VC90")
ELSE(CMAKE_GENERATOR MATCHES "Visual Studio 9")
  IF(CMAKE_GENERATOR MATCHES "Visual Studio 8")
    SET(VS_VERSION "VC80")
  ELSE(CMAKE_GENERATOR MATCHES "Visual Studio 8")
    MESSAGE( STATUS "Automatically generating path to redistributable libraries is not supported for ${CMAKE_GENERATOR}.  Please edit InnoScript (*.iss) files in ${SCIRUN_CORE_BINARY_DIR}/InstallUtils/windows after cmake generates build files.")
  ENDIF(CMAKE_GENERATOR MATCHES "Visual Studio 8")
ENDIF(CMAKE_GENERATOR MATCHES "Visual Studio 9")

IF(CMAKE_SIZEOF_VOID_P MATCHES 8)
  SET(WIN_REDIST "amd64")
  SET(SETUP_ARCH "x64")
  SET(WIN_SCIRUN_BUILD_TYPE "win64")
ELSE(CMAKE_SIZEOF_VOID_P MATCHES 8)
  SET(WIN_REDIST "x86")
  SET(SETUP_ARCH "") # if blank, setup runs in 32-bit mode
  SET(WIN_SCIRUN_BUILD_TYPE "win32")
ENDIF(CMAKE_SIZEOF_VOID_P MATCHES 8)

FOREACH(BUILD_TYPE Debug Release RelWithDebInfo MinSizeRel) 
  IF(BUILD_BIOMESH3D)
    SET(DISABLE_BIOMESH3D_SCRIPTS "")
  ELSE(BUILD_BIOMESH3D)
    SET(DISABLE_BIOMESH3D_SCRIPTS ";")
  ENDIF(BUILD_BIOMESH3D)

  IF(${ITK_BUILD_SHARED_LIBS} MATCHES "SHARED")
    SET(DISABLE_ITK "")
  ELSE(${ITK_BUILD_SHARED_LIBS} MATCHES "SHARED")
    SET(DISABLE_ITK ";")
  ENDIF(${ITK_BUILD_SHARED_LIBS} MATCHES "SHARED")

  IF(BUILD_TYPE MATCHES "Debug")
    SET(DEBUG_BUILD "Debug")
    SET(DEBUG_REDIST "Debug_NonRedist")
  ELSE(BUILD_TYPE MATCHES "Debug")
    SET(DEBUG_BUILD "") 
    SET(DEBUG_REDIST "")
    SET(DISABLE_PDB ";")
  ENDIF(BUILD_TYPE MATCHES "Debug")

  SET(CRT_REDIST "Microsoft.${VS_VERSION}.${DEBUG_BUILD}CRT")
  SET(MFC_REDIST "Microsoft.${VS_VERSION}.${DEBUG_BUILD}MFC")
  SET(PATH_TO_CRT_REDIST "${PATH_TO_VC}/redist/${DEBUG_REDIST}/${WIN_REDIST}/${CRT_REDIST}")
  SET(PATH_TO_MFC_REDIST "${PATH_TO_VC}/redist/${DEBUG_REDIST}/${WIN_REDIST}/${MFC_REDIST}")
  FILE(TO_NATIVE_PATH "${PATH_TO_CRT_REDIST}" PATH_TO_CRT_REDIST_NATIVE)
  FILE(TO_NATIVE_PATH "${PATH_TO_MFC_REDIST}" PATH_TO_MFC_REDIST_NATIVE)
  FILE(TO_NATIVE_PATH "${ITK_DIR}" ITK_DIR_NATIVE)

  CONFIGURE_FILE (
    "${SCIRUN_CORE_SOURCE_DIR}/InstallUtils/windows/scirun.iss.in"
    "${SCIRUN_CORE_BINARY_DIR}/InstallUtils/windows/scirun_${BUILD_TYPE}.iss"
  )

  IF(BUILD_SEG3D)
    IF(BUILD_TYPE MATCHES "Debug")
      SET(SEG3D_BASE_FILENAME "Seg3DDebug_${SEG3D_VERSION_MAJOR}.  ${SEG3D_VERSION_MINOR}")
    ELSE(BUILD_TYPE MATCHES "Debug")
      SET(SEG3D_BASE_FILENAME "Seg3D_${SEG3D_VERSION_MAJOR}.${SEG3D_VERSION_MINOR}")
      SET(DISABLE_WX_PDB ";")
    ENDIF(BUILD_TYPE MATCHES "Debug")

    IF(WITH_WXWIDGETS)
      FILE(TO_NATIVE_PATH "${WX_LIB_DIR}" WX_LIB_DIR_NATIVE)
    ELSE(WITH_WXWIDGETS)
      SET(DISABLE_WX_LIBS ";")
      SET(DISABLE_WX_PDB ";")
    ENDIF(WITH_WXWIDGETS)

    CONFIGURE_FILE (
      "${SCIRUN_CORE_SOURCE_DIR}/InstallUtils/windows/seg3d.iss.in"
      "${SCIRUN_CORE_BINARY_DIR}/InstallUtils/windows/seg3d_${BUILD_TYPE}.iss"
    )
  ENDIF(BUILD_SEG3D)
ENDFOREACH(BUILD_TYPE Debug Release RelWithDebInfo MinSizeRel) 

CONFIGURE_FILE (
  "${SCIRUN_CORE_SOURCE_DIR}/InstallUtils/windows/buildReleaseInstaller.bat.in"
  "${SCIRUN_CORE_BINARY_DIR}/buildReleaseInstaller.bat"
)
