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

# 
SET(IWIDGETS_VERSION "4.0.1" CACHE INTERNAL "IWidgets version")


#################################################################
# Building directories with tcl files

SET(ITCL_WIDGETS_DIRECTORY
      ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Iwidgets${IWIDGETS_VERSION}
    CACHE INTERNAL
    "Path to iwidgets scripts")

FILE(MAKE_DIRECTORY ${ITCL_WIDGETS_DIRECTORY})

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/generic/tclIndex ${ITCL_WIDGETS_DIRECTORY}/tclIndex COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/iwidgets.tcl ${ITCL_WIDGETS_DIRECTORY}/iwidgets.tcl COPYONLY)

SET(VERSION "${IWIDGETS_VERSION}")
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/pkgIndex.tcl.in ${ITCL_WIDGETS_DIRECTORY}/pkgIndex.tcl)

FILE(GLOB itk_files "${CMAKE_CURRENT_SOURCE_DIR}/generic/*.itk")
FOREACH(itk_file ${itk_files})
  GET_FILENAME_COMPONENT(itk_name ${itk_file} NAME)
  CONFIGURE_FILE(${itk_file}
                 ${ITCL_WIDGETS_DIRECTORY}/${itk_name}
                 COPYONLY)
ENDFOREACH(itk_file ${tik_files})

FILE(GLOB itcl_files "${CMAKE_CURRENT_SOURCE_DIR}/generic/*.itcl")
FOREACH(itcl_file ${itcl_files})
  GET_FILENAME_COMPONENT(itcl_name ${itcl_file} NAME)
  CONFIGURE_FILE(${itcl_file}
                 ${ITCL_WIDGETS_DIRECTORY}/${itcl_name}
                 COPYONLY)
ENDFOREACH(itcl_file ${itcl_files})
