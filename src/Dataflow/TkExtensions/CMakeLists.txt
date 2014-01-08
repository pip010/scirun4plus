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

# CMakeLists.txt for Dataflow/TkExtensions

IF(WITH_X11)
      ADD_DEFINITIONS(-DGLEW_APPLE_GLX)
ENDIF(WITH_X11)


SET(Dataflow_TkExtensions_SRCS
  tk3d2.c
  tkBevel.c
  tkCanvBLine.c
  tkCursor.c
  tkOpenGL.c
  tk3daux.c
  tkRange.c
  tkUnixRange.c
)

SCIRUN_ADD_LIBRARY(Dataflow_TkExtensions ${Dataflow_TkExtensions_SRCS})

TARGET_LINK_LIBRARIES(Dataflow_TkExtensions
  Core_Thread
  Core_Geom
  ${SCI_BLT_LIBRARY}
  ${SCI_ITK_LIBRARY}
  ${SCI_ITCL_LIBRARY}
  ${SCI_TK_LIBRARY}
  ${SCI_TCL_LIBRARY}
  ${GL_LIBRARY}
  ${SCI_GLEW_LIBRARY}
  ${M_LIBRARY}
  ${X_LIBRARY}
)

IF(APPLE_LEOPARD)
  SET_TARGET_PROPERTIES(Dataflow_TkExtensions PROPERTIES LINK_FLAGS  "${GL_LINK} ${X_LINK}")
ENDIF(APPLE_LEOPARD)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Dataflow_TkExtensions)
ENDIF(BUILD_SHARED_LIBS)

FOREACH(definition ${SCI_TK_DEFINITIONS})
  ADD_DEFINITIONS("-D${definition}")
ENDFOREACH(definition ${SCI_TK_DEFINITIONS})

IF(WIN32)
  SET_TARGET_PROPERTIES(Dataflow_TkExtensions
                        PROPERTIES
                          COMPILE_FLAGS "-fp:strict -Dinline=__inline -DBUILD_tk -DBUILD_ttk -DSUPPORT_CONFIG_EMBEDDED -D_MBCS"
                          DEFINE_SYMBOL _USRDLL
                          DEFINE_SYMBOL _WINDLL
                    	)
ENDIF(WIN32)


