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

# This is the root SCIRun CMakeLists file.
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(SCIRUN_DATAFLOW)

# TODO: temporary -- test on *nix
#ADD_DEFINITIONS(-DGLEW_STATIC)
#ADD_DEFINITIONS(-DSLIVR_STATIC)

SET(SCIRUN_TCL_VERSION_MAJOR "1" CACHE INTERNAL "SCIRun TCL package major version number." FORCE)
SET(SCIRUN_TCL_VERSION_MINOR "0" CACHE INTERNAL "SCIRun TCL package minor version number." FORCE)
SET(SCIRUN_TCL_VERSION "${SCIRUN_TCL_VERSION_MAJOR}.${SCIRUN_TCL_VERSION_MINOR}" CACHE INTERNAL "SCIRun TCL package version string." FORCE)

#-----------------------------------------------------------------------------
# Test for some required system information.
INCLUDE (CMakeBackwardCompatibilityC)

#-----------------------------------------------------------------------------
# Configure testing support.  This should be done before any
# MESSAGE(FATAL_ERROR ...) commands are invoked.
INCLUDE( CTest )
IF(BUILD_TESTING)
  ENABLE_TESTING()
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)

#----------------------------------------------------------------------------
# Output directories.
IF(NOT BUILD_FROM_TOP_LEVEL)
  SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${SCIRUN_DATAFLOW_BINARY_DIR}/lib CACHE INTERNAL "Single output directory for building all libraries.")
  SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SCIRUN_DATAFLOW_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
  #output directory for installing all include files.")
  MARK_AS_ADVANCED(CMAKE_LIBRARY_OUTPUT_DIRECTORY CMAKE_RUNTIME_OUTPUT_DIRECTORY)
ENDIF(NOT BUILD_FROM_TOP_LEVEL)

# DV Fix when fixing separate building option
IF(WITH_ITK)
  IF(NOT ITK_DIR)
    SET(ITK_DIR ${SCIRUN_ITK_DIR})
  ENDIF(NOT ITK_DIR)

  FIND_PACKAGE( ITK )
  IF(ITK_FOUND)
    # THIS HACK SHOULD PREVENT CMAKE FROM SPILLING ALL THE ITK SETTINGS
    # INTO THE SCIRUN CORE SETTINGS
    SET(ITK_BUILD_SETTINGS_FILE "") 

    IF(ITK_USE_FILE)
      INCLUDE(${ITK_USE_FILE})
    ELSE(ITK_USE_FILE)
      INCLUDE(${SCI_ITK_USE_FILE})
    ENDIF(ITK_USE_FILE)
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR "SCIRun Core build with the Insight Toolkit. Please set ITK_DIR.")
  ENDIF(ITK_FOUND)
ENDIF(WITH_ITK)

SET(scirun_obj_dataflow_dir ${SCIRUN_DATAFLOW_BINARY_DIR})

SET(IS_OSX "no")
SET(IS_WIN "no")
SET(IS_AIX "no")

# Add correct flags for win32
IF(WIN32)
  IF(NOT BUILD_SHARED_LIBS)
    ADD_DEFINITIONS(-DBUILD_SCIRUN_STATIC)
  ENDIF(NOT BUILD_SHARED_LIBS)

#  This defines M_PI for windows
  ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
  ADD_DEFINITIONS(-D_WIN32_WINNT=0x0500)
  # win32 l3 warnings are a little ridiculous
  SET(CMAKE_CXX_WARNING_LEVEL 1)
  IF(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    STRING(REGEX REPLACE "/W[0-4]" "/W1"
            CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  ELSE(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W1")
  ENDIF(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
ENDIF(WIN32)

ADD_SUBDIRECTORY(Comm)
ADD_SUBDIRECTORY(Constraints)
ADD_SUBDIRECTORY(GuiInterface)
ADD_SUBDIRECTORY(Modules)
ADD_SUBDIRECTORY(Network)
ADD_SUBDIRECTORY(TkExtensions)
ADD_SUBDIRECTORY(Widgets)

###############################################################################
# BUILD PACKAGE LIST

STRING(REGEX REPLACE "[,]" ";" PACKAGE_LIST ${LOAD_PACKAGE})
SET(INCLUDED_PACKAGE_INCLUDES FALSE)

FOREACH( pkg ${PACKAGE_LIST})
  # only do anything for packages other than SCIRun
  IF(NOT ${pkg} MATCHES "SCIRun" AND NOT HANDLED_PACKAGE_${pkg})

    IF(NOT INCLUDED_PACKAGE_INCLUDES)
      INCLUDE_DIRECTORIES(${SCIRUN_CORE_SOURCE_DIR}/Packages)
      SET(INCLUDED_PACKAGE_INCLUDES TRUE)
    ENDIF(NOT INCLUDED_PACKAGE_INCLUDES)

    # Add the package subdirectory in CMake
    ADD_SUBDIRECTORY(${SCIRUN_CORE_SOURCE_DIR}/Packages/${pkg} ${pkg})

    # Install the relevant source directories
    SCIRUN_INSTALL_PACKAGE(${pkg})
    
  ENDIF(NOT ${pkg} MATCHES "SCIRun" AND NOT HANDLED_PACKAGE_${pkg})
ENDFOREACH( pkg )


###############################################################################
# BUILDING TCL INDEX FILES

SET(SCIRUN_TCL_DIRECTORY
      ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/sciruntcl${SCIRUN_TCL_VERSION}
    CACHE INTERNAL
      "Path to SCIRun tcl scripts"
    FORCE)

FILE(MAKE_DIRECTORY ${SCIRUN_TCL_DIRECTORY})

# Copy SCIRun Dataflow GUI tcl files from source to binary directory
FILE(GLOB tcl_files "${SCIRUN_CORE_SOURCE_DIR}/Dataflow/GUI/*.tcl")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/GUI/${tcl_name}
                 ${SCIRUN_TCL_DIRECTORY}/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

SET(CREATE_TCL_PACKAGE ${SCIRUN_TCL_DIRECTORY})

FOREACH( pkg ${PACKAGE_LIST})
  IF(NOT ${pkg} MATCHES "SCIRun" AND NOT HANDLED_PACKAGE_${pkg})

  #SET(CREATE_TCL_PACKAGE ${SCIRUN_CORE_SOURCE_DIR}/Packages/${pkg}/Dataflow/GUI)

  FILE(GLOB tcl_files "${SCIRUN_CORE_SOURCE_DIR}/Packages/${pkg}/Dataflow/GUI/*.tcl")
  FOREACH(tcl_file ${tcl_files})
    GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
    CONFIGURE_FILE(${tcl_file}
                 ${SCIRUN_TCL_DIRECTORY}/${tcl_name}
                 COPYONLY)
  ENDFOREACH(tcl_file ${tcl_files})

  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/GUI/PowerAppBase.tcl
                             ${SCIRUN_TCL_DIRECTORY}/PowerAppBase.tcl
                             COPYONLY)
 
  ENDIF(NOT ${pkg} MATCHES "SCIRun" AND NOT HANDLED_PACKAGE_${pkg})
ENDFOREACH( pkg )

CONFIGURE_FILE (
  "${CMAKE_CURRENT_SOURCE_DIR}/GUI/cmake_createTclIndex.tcl.in"
  "${SCIRUN_TCL_DIRECTORY}/createSCIRunTclIndex.tcl"
)

GET_TARGET_PROPERTY(SCI_TCLSH_LOCATION ${SCI_TCLSH} LOCATION)
ADD_CUSTOM_TARGET(create_tclIndex_SCIRun ALL
  COMMAND ${SCI_TCLSH_LOCATION} ${SCIRUN_TCL_DIRECTORY}/createSCIRunTclIndex.tcl
  DEPENDS ${SCI_TCLSH}
)  

###############################################################################

# Create powerapps
STRING (REGEX MATCH "Teem" WITH_TEEM_PACKAGE ${LOAD_PACKAGE})
STRING (REGEX MATCH "BioPSE" WITH_BIOPSE_PACKAGE ${LOAD_PACKAGE})

IF(WITH_TEEM_PACKAGE MATCHES "Teem")
  # BioTensor
  IF(WIN32)
    FOREACH(BUILD_TYPE Debug Release)
      CONFIGURE_FILE(${SCIRUN_CORE_SOURCE_DIR}/Packages/Teem/nets/PowerApps/cmake_winBioTensor
                                 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${BUILD_TYPE}/BioTensor.bat)
    ENDFOREACH(BUILD_TYPE Debug Release)
  ELSE(WIN32)
    CONFIGURE_FILE(${SCIRUN_CORE_SOURCE_DIR}/Packages/Teem/nets/PowerApps/BioTensor.app
                               ${SCIRUN_TCL_DIRECTORY}/BioTensor.app
                               COPYONLY)
    CONFIGURE_FILE(${SCIRUN_CORE_SOURCE_DIR}/Packages/Teem/nets/PowerApps/cmake_BioTensor
                               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BioTensor)

    ADD_CUSTOM_TARGET(create_BioTensor
                      ALL
                      chmod +x ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BioTensor)
  ENDIF(WIN32)
ENDIF(WITH_TEEM_PACKAGE MATCHES "Teem")
