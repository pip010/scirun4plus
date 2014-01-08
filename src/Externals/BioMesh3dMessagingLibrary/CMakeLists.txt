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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(BioMesh3dMessage)
###########################################
# Setup variables that other parts of the
# make system will use to link to this
# library.
###########################################

SET(BIOMESH3D_MESSAGE_LIBRARY BioMesh3dMessage CACHE INTERNAL "MediaMessage library name." FORCE)
SET(BIOMESH3D_MESSAGE_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "MediaMessage include." FORCE)

###########################################
# Source files
###########################################

set( BioMesh3dMessage_SRCS
     ClientMessageComposer.cpp
     ClientMessageHandler.cpp
     SCIRunMessageComposer.cpp
     SCIRunMessageHandler.cpp
     ServerMessageComposer.cpp
     ServerMessageHandler.cpp
     BrokerMessageComposer.cpp
     BrokerMessageHandler.cpp
)

set( BioMesh3dMessage_HDRS
     ClientMessageComposer.h
     ClientMessageHandler.h
     SCIRunMessageComposer.h
     SCIRunMessageHandler.h
     ServerMessageComposer.h
     ServerMessageHandler.h
     BrokerMessageComposer.h
     BrokerMessageHandler.h
)

###########################################
# Create the library
###########################################

INCLUDE_DIRECTORIES(${SCI_BOOST_INCLUDE} ${SCI_CRYPTO_INCLUDE} ${SCI_MEDIAMESSAGE_INCLUDE})
ADD_LIBRARY(${BIOMESH3D_MESSAGE_LIBRARY} STATIC ${BioMesh3dMessage_SRCS} ${BioMesh3dMessage_HDRS})

TARGET_LINK_LIBRARIES( BioMesh3dMessage ${SCI_BOOST_LIBRARY} ${SCI_CRYPTO_LIBRARY} 
	${MEDIA_MESSAGE_LIBRARY} ${MEDIA_COMMUNICATION_LIBRARY} MML_Core_Utils MML_Core_Codec) 

###########################################
# Linux needs the PIC flag
###########################################

IF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)
  SET_TARGET_PROPERTIES( ${BIOMESH3D_MESSAGE_LIBRARY} PROPERTIES COMPILE_FLAGS "-fPIC" )
ENDIF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)
