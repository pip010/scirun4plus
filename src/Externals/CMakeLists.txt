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

MESSAGE(STATUS "Configuring External Packages")

##################################################
## Add zlib

MESSAGE(STATUS "Configuring zlib")

ADD_SUBDIRECTORY(zlib)

##################################################
## Add libpng

MESSAGE(STATUS "Configuring libpng")

ADD_SUBDIRECTORY(libpng)

##################################################
## Add Teem

## Teem depends on libpng and zlib, so make sure that their projects are
## included before Teem.

MESSAGE(STATUS "Configuring teem")

ADD_SUBDIRECTORY(teem)

##################################################
## Add glew

MESSAGE(STATUS "Configuring glew")

ADD_SUBDIRECTORY(glew)

##################################################
## Add SLIVR

SET(SCIRUN_SLIVR_BUILD ON CACHE INTERNAL "Build slivr as a sub-project.")

## slivr depends on teem and glew, so make sure that their projects are
## included before slivr.

MESSAGE(STATUS "Configuring slivr")

ADD_SUBDIRECTORY(slivr)

##################################################
## Add FreeType

MESSAGE(STATUS "Configuring freetype")

ADD_SUBDIRECTORY(freetype)

##################################################
## Add libXML2

MESSAGE(STATUS "Configuring libxml2")

ADD_SUBDIRECTORY(libxml2)

##################################################
# Add TCL/Tk/iTCL/iTK/iWidgets

IF(BUILD_DATAFLOW)
  ADD_SUBDIRECTORY(tcltk)
ENDIF(BUILD_DATAFLOW)

##################################################
# Add tetgen

IF(WITH_TETGEN OR BUILD_BIOMESH3D)
  MESSAGE(STATUS "Configuring tetgen")
  ADD_SUBDIRECTORY(tetgen)
ENDIF(WITH_TETGEN OR BUILD_BIOMESH3D)

##################################################
# Add particle-system-mm and vispack

IF(BUILD_BIOMESH3D)
  MESSAGE(STATUS "Configuring particle system")
  ADD_SUBDIRECTORY(particle-system-mm)

  MESSAGE(STATUS "Configuring vispack")
  ADD_SUBDIRECTORY(vispack)
ENDIF(BUILD_BIOMESH3D)

##################################################
# Add ITK

IF(WITH_ITK OR BUILD_SEG3D)
  MESSAGE(STATUS "Configuring itk")
  ADD_SUBDIRECTORY(itk)
ENDIF(WITH_ITK OR BUILD_SEG3D)

##################################################
# Add Boost

MESSAGE(STATUS "Configuring Boost")
ADD_SUBDIRECTORY(boost)

IF(BUILD_BIOMESH3D_REMOTE_SUPPORT)

#################################################
# Add FFMPEG Library

ADD_SUBDIRECTORY(ffmpeg)

#################################################
# Add openssl

ADD_SUBDIRECTORY(ssl)

#################################################
# Add libvpx
ADD_SUBDIRECTORY(libvpx)

#################################################
# Add MediaMessagingLibrary

ADD_SUBDIRECTORY(MediaMessagingLibrary/src)

#################################################
# Add MediaMessagingLibrary

INCLUDE_DIRECTORIES(${MEDIA_MESSAGE_INCLUDE})

ADD_SUBDIRECTORY(BioMesh3dMessagingLibrary)

ENDIF(BUILD_BIOMESH3D_REMOTE_SUPPORT)

##################################################
# Add Netlib

IF(BUILD_NETLIB)
  MESSAGE(STATUS "Configuring Netlib")
  ADD_SUBDIRECTORY(netlib)
ENDIF(BUILD_NETLIB)

##################################################
# Add Google Test

IF(BUILD_TESTING)
  MESSAGE(STATUS "Configuring Google Test")
  ADD_SUBDIRECTORY(googletest)
ENDIF(BUILD_TESTING)
