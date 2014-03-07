#   zlib source code license
#
#   Copyright (C) 1995-2005 Jean-loup Gailly and Mark Adler
# 
#   This software is provided 'as-is', without any express or implied
#   warranty.  In no event will the authors be held liable for any damages
#   arising from the use of this software.
# 
#   Permission is granted to anyone to use this software for any purpose,
#   including commercial applications, and to alter it and redistribute it
#   freely, subject to the following restrictions:
# 
#   1. The origin of this software must not be misrepresented; you must not
#      claim that you wrote the original software. If you use this software
#      in a product, an acknowledgment in the product documentation would be
#      appreciated but is not required.
#   2. Altered source versions must be plainly marked as such, and must not be
#      misrepresented as being the original software.
#   3. This notice may not be removed or altered from any source distribution.
# 
#   Jean-loup Gailly        Mark Adler
#   jloup@gzip.org          madler@alumni.caltech.edu
# 
# 
#   The data format used by the zlib library is described by RFCs (Request for
#   Comments) 1950 to 1952 in the files http://www.ietf.org/rfc/rfc1950.txt
#   (zlib format), rfc1951.txt (deflate format) and rfc1952.txt (gzip format).

# CMake build scripts provided by Scientific Computing and Imaging Institute
# Copyright (c) 2009 Scientific Computing and Imaging Institute, University of Utah.


CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(VERSION 2.6)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(ZLIB)

IF(BUILD_STATIC_THIRDPARTY)
  IF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)
    SET(ENABLE_PIC ON)
  ENDIF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)
  SET(BUILD_SHARED_LIBS OFF)
ENDIF(BUILD_STATIC_THIRDPARTY)

#-----------------------------------------------------------------------------
# Test for some required system information.
INCLUDE (CMakeBackwardCompatibilityC)

# source files for zlib
SET(ZLIB_SRCS
      adler32.c
      compress.c
      crc32.c
      deflate.c
      gzio.c
      infback.c
      inffast.c
      inflate.c
      inftrees.c
      trees.c
      uncompr.c
      zutil.c
    )

# for windows add the .def and .rc files to the source list
# if building shared libs
IF(WIN32)
  IF(MSVC)
    # Bind application to the current version of the VC libraries
    ADD_DEFINITIONS(-D_BIND_TO_CURRENT_VCLIBS_VERSION) 

    # MSVC level 3 warnings are a little ridiculous
    SET(CMAKE_C_WARNING_LEVEL 1)
    IF(CMAKE_C_FLAGS MATCHES "/W[0-4]")
      STRING(REGEX REPLACE "/W[0-4]" "/W1"
              CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    ELSE(CMAKE_C_FLAGS MATCHES "/W[0-4]")
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W1")
    ENDIF(CMAKE_C_FLAGS MATCHES "/W[0-4]")
  ENDIF(MSVC)

  IF(BUILD_SHARED_LIBS)
    IF(NOT UNIX AND NOT BORLAND AND NOT MINGW)
      SET(ZLIB_SRCS ${ZLIB_SRCS} win32/zlib.def win32/zlib1.rc)
    ENDIF(NOT UNIX AND NOT BORLAND AND NOT MINGW)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(WIN32)

IF(WIN32)
  SET(SCI_ZLIB_LIBRARY zlib1 CACHE INTERNAL "Zlib library name." FORCE)
ELSE(WIN32)
  SET(SCI_ZLIB_LIBRARY zlib CACHE INTERNAL "Zlib library name." FORCE)
ENDIF(WIN32)
SET(SCI_ZLIB_INCLUDE
      ${CMAKE_CURRENT_SOURCE_DIR}
    CACHE INTERNAL
      "Zlib include directories."
    FORCE)

INCLUDE_DIRECTORIES(${SCI_ZLIB_INCLUDE})
ADD_LIBRARY(${SCI_ZLIB_LIBRARY} ${ZLIB_SRCS})

IF(WIN32)
  IF(NOT UNIX AND NOT BORLAND AND NOT MINGW)
    SET_TARGET_PROPERTIES(${SCI_ZLIB_LIBRARY}
                          PROPERTIES
                            COMPILE_FLAGS "/Ob1 /Gy /GF"
                            LINK_FLAGS "/DYNAMICBASE:NO")
  ENDIF(NOT UNIX AND NOT BORLAND AND NOT MINGW)
  IF(BUILD_SHARED_LIBS)
    SET_TARGET_PROPERTIES(${SCI_ZLIB_LIBRARY}
                          PROPERTIES
                            DEFINE_SYMBOL
                            ZLIB_DLL)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(WIN32)

IF(ENABLE_PIC)
  SET_TARGET_PROPERTIES(${SCI_ZLIB_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fPIC")
ENDIF(ENABLE_PIC)


