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
# CMake build scripts provided by Scientific Computing and Imaging Institute
# Copyright (c) 2009 Scientific Computing and Imaging Institute, University of Utah.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(VERSION 2.6)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(tcltk)

# SOME HOW CMAKE DOES NOT GET IT RIGHT
IF(APPLE)
IF(CMAKE_OSX_ARCHITECTURES MATCHES .*64)
  SET(CMAKE_SIZEOF_VOID_P 8)
ENDIF(CMAKE_OSX_ARCHITECTURES MATCHES .*64)
ENDIF(APPLE)

IF(BUILD_STATIC_THIRDPARTY)
  IF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)
    SET(ENABLE_PIC ON)
  ENDIF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)
  SET(BUILD_SHARED_LIBS OFF)
ENDIF(BUILD_STATIC_THIRDPARTY)

#-----------------------------------------------------------------------------
# Test for some required system information.
INCLUDE (CMakeBackwardCompatibilityC)
INCLUDE(CheckCSourceCompiles)
INCLUDE(CheckFunctionExists)

# ADD TCL TO PROJECT
MESSAGE(STATUS "Configuring TCL")

ADD_SUBDIRECTORY(tcl8.6b1)

MESSAGE(STATUS "Configuring TK")

ADD_SUBDIRECTORY(tk8.6b1)

MESSAGE(STATUS "Configuring ITCL")

ADD_SUBDIRECTORY(itcl4.0b3)

MESSAGE(STATUS "Configuring ITK")

ADD_SUBDIRECTORY(itk4.0b3)

MESSAGE(STATUS "Configuring iWidgets")

ADD_SUBDIRECTORY(iwidgets4.0.1)

MESSAGE(STATUS "Configuring BLT")

ADD_SUBDIRECTORY(blt2.4z)
