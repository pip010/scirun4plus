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

# CMakeLists.txt for Dataflow/GuiInterface

SET(Dataflow_GuiInterface_SRCS
  GuiCallback.cc
  GuiContext.cc
  GuiGeom.cc
  GuiVar.cc
  GuiView.cc
  TCLInterface.cc
  TCLThread.cc
  UIvar.cc
  TkOpenGLContext.cc
)

SCIRUN_ADD_LIBRARY(Dataflow_GuiInterface ${Dataflow_GuiInterface_SRCS})

TARGET_LINK_LIBRARIES(Dataflow_GuiInterface
  Core_Exceptions
  Core_Util
  Core_Thread
  Dataflow_TkExtensions
  Core_Math
  Core_Geom
  Core_Datatypes
  ${ITK_LIBRARY}
  ${ITCL_LIBRARY}
  ${TCL_LIBRARY}
  ${X_LIBRARY}
  ${SCI_TEEM_LIBRARY}
  ${GL_LIBRARY}
)

IF(APPLE_LEOPARD)
  SET_TARGET_PROPERTIES(Dataflow_GuiInterface PROPERTIES LINK_FLAGS  "${GL_LINK} ${X_LINK}")
ENDIF(APPLE_LEOPARD)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Dataflow_GuiInterface)
ENDIF(BUILD_SHARED_LIBS)

FOREACH(definition ${SCI_TK_DEFINITIONS})
  ADD_DEFINITIONS("-D${definition}")
ENDFOREACH(definition ${SCI_TK_DEFINITIONS})

IF(WIN32)
  SET_TARGET_PROPERTIES(Dataflow_GuiInterface
                        PROPERTIES
                          COMPILE_FLAGS "-fp:strict -Dinline=__inline -DTCL_THREADS=1 -DSUPPORT_CONFIG_EMBEDDED -D_MBCS"
                          DEFINE_SYMBOL _USRDLL
                          DEFINE_SYMBOL _WINDLL
                    	)
ENDIF(WIN32)

