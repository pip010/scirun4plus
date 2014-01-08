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

# This is the root ITK CMakeLists file.  
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(VERSION 2.6)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(BIOMESH3D)

MESSAGE(STATUS "Configuring BioMesh3D")

SET(BIOMESH_VERSION_MAJOR 0  CACHE INTERNAL "BioMesh3D major version number." FORCE)
SET(BIOMESH_VERSION_MINOR 1  CACHE INTERNAL "BioMesh3D minor version number." FORCE)
SET(BIOMESH_VERSION_PATCH 0  CACHE INTERNAL "BioMesh3D patch number." FORCE)
SET(BIOMESH_VERSION_STRING "${BIOMESH_VERSION_MAJOR}.${BIOMESH_VERSION_MINOR}"  CACHE INTERNAL "BioMesh3D version string." FORCE)

#FIND_PACKAGE(PythonInterp)
#
#IF(NOT PYTHONINTERP_FOUND)
#  MESSAGE(FATAL_ERROR "Python is required to run the BioMesh3D pipeline scripts.")
#ENDIF(NOT PYTHONINTERP_FOUND)


## TODO: we also need this for install tools, so should be moved to top-level cmake file
STRING(REPLACE "/src" "" SCIRUN_TOP_DIR ${SCIRUN_CORE_SOURCE_DIR})
IF(WIN32)
  FILE(TO_NATIVE_PATH "${SCIRUN_CORE_BINARY_DIR}" FEMESHER_SCIRUN_BIN)
  FILE(TO_NATIVE_PATH "${SCIRUN_TOP_DIR}" FEMESHER_SCIRUN_HOME)
ELSE(WIN32)
  SET(FEMESHER_SCIRUN_BIN ${SCIRUN_CORE_BINARY_DIR})
  SET(FEMESHER_SCIRUN_HOME ${SCIRUN_TOP_DIR})
ENDIF(WIN32)

# Python script directory
FILE(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FEMesher)

# Copy Python scripts

FILE(GLOB python_files "${BIOMESH3D_SOURCE_DIR}/scripts/*.py")

FOREACH(python_file ${python_files})
  GET_FILENAME_COMPONENT(python_name ${python_file} NAME)
  CONFIGURE_FILE(${python_file}
                 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FEMesher/${python_name}
                 COPYONLY)
ENDFOREACH(python_file ${python_files})

# Copy SCIRun networks

FILE(GLOB srn_files "${SCIRUN_CORE_SOURCE_DIR}/nets/MeshingPipeline/*.srn")

FOREACH(srn_file ${srn_files})
  GET_FILENAME_COMPONENT(srn_name ${srn_file} NAME)
  CONFIGURE_FILE(${srn_file}
                 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/FEMesher/${srn_name}
                 COPYONLY)
ENDFOREACH(srn_file ${srn_files})

