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
INCLUDE(CheckSymbolExists)

SET(BLT_VERSION "2.4" CACHE INTERNAL "BLT version")
SET(BLT_MAJOR_VERSION "2")
SET(BLT_MINOR_VERSION "4")
SET(BLT_PATCH_LEVEL "")

#################################################################
# Setting up include files and source files

SET(BLT_GRAPH_SRCS
    src/bltGrAxis.c
		src/bltGrBar.c
		src/bltGrElem.c
		src/bltGrGrid.c
		src/bltGrHairs.c
		src/bltGrLegd.c
		src/bltGrLine.c
		src/bltGrMarker.c
		src/bltGrMisc.c
		src/bltGrPen.c
		src/bltGrPs.c
		src/bltGraph.c   
    )

SET(BLT_TCL_ONLY_SRCS
    src/bltAlloc.c
		src/bltArrayObj.c
		src/bltBgexec.c
		src/bltChain.c
		src/bltDebug.c
		src/bltHash.c
		src/bltList.c
		src/bltNsUtil.c
		src/bltParse.c
		src/bltPool.c
		src/bltSpline.c
		src/bltSwitch.c
		src/bltTree.c
		src/bltTreeCmd.c
		src/bltUtil.c
		src/bltVecCmd.c
		src/bltVecMath.c
		src/bltVecObjCmd.c
		src/bltVector.c
		src/bltWatch.c
   )

SET(BLT_GENERIC_SRCS
    src/bltTed.c
    src/tkButton.c
    src/tkFrame.c
    src/bltScrollbar.c
		src/bltBeep.c
		src/bltBind.c
		src/bltBitmap.c
		src/bltBusy.c
		src/bltCanvEps.c
		src/bltConfig.c
		src/bltContainer.c
		src/bltDragdrop.c
		src/bltHierbox.c
		src/bltHtext.c
		src/bltImage.c
		src/bltPs.c
		src/bltTable.c
		src/bltTabnotebook.c
		src/bltTabset.c
		src/bltText.c
		src/bltTile.c
		src/bltTreeView.c
		src/bltTreeViewCmd.c
		src/bltTreeViewEdit.c
		src/bltTreeViewColumn.c
		src/bltTreeViewStyle.c
    src/bltWindow.c
		src/bltObjConfig.c
		src/bltWinop.c
    src/bltInit.c
  )

SET(BLT_UNIX_SRCS
    src/bltColor.c
    src/bltCutbuffer.c
		src/bltUnixPipe.c
    src/bltUnixImage.c
    src/bltUnixDnd.c
  )

SET(BLT_WIN_SRCS
		src/bltWinDraw.c
		src/bltWinPrnt.c
		src/bltWinImage.c
		src/bltWinPipe.c
		src/bltWinUtil.c
		src/bltWinDde.c
		src/pure_api.c
   )

IF (UNIX)
  # X11 version (so we can run 64bit on OSX as well)
  SET(BLT_SRCS ${BLT_GRAPH_SRCS} ${BLT_TCL_ONLY_SRCS} ${BLT_GENERIC_SRCS} ${BLT_UNIX_SRCS})
  FIND_PACKAGE( X11 )
  INCLUDE_DIRECTORIES(
    src
    ${SCI_TCL_INCLUDE}
    ${SCI_TK_INCLUDE}
    ${X11_INCLUDE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
  )
ELSE(UNIX)
  # Windows version

  IF(WIN32)
    SET(BLT_SRCS ${BLT_GRAPH_SRCS} ${BLT_TCL_ONLY_SRCS} ${BLT_GENERIC_SRCS} ${BLT_WIN_SRCS})

    IF(BUILD_SHARED_LIBS)
      # Default machine type is IX86.
      # MACHINE=(IX86|IA64|AMD64|ALPHA)
      IF(CMAKE_SIZEOF_VOID_P MATCHES 8)
        SET(MACHINE "AMD64")
      ELSE(CMAKE_SIZEOF_VOID_P MATCHES 8)
        SET(MACHINE "IX86")
      ENDIF(CMAKE_SIZEOF_VOID_P MATCHES 8)
    ENDIF(BUILD_SHARED_LIBS)

    INCLUDE_DIRECTORIES(
      src
      ${SCI_TCL_INCLUDE}
      ${SCI_TK_INCLUDE}
      ${CMAKE_CURRENT_BINARY_DIR}
    )

  ENDIF(WIN32)
ENDIF(UNIX)

#################################################################
# Setting up definitions

ADD_DEFINITIONS(-DPACKAGE_NAME="BLT")
ADD_DEFINITIONS(-DPACKAGE_TARNAME="BLT")
ADD_DEFINITIONS(-DPACKAGE_VERSION="${BLT_VERSION}")
ADD_DEFINITIONS(-DPACKAGE_STRING="BLT ${BLT_VERSION}")
ADD_DEFINITIONS(-DPACKAGE_BUGREPORT="")

CHECK_C_SOURCE_COMPILES(
                         "
                           #include <stddef.h>
                           #include <stdlib.h>

                           int main()
                           {
                             return 0;
                           }
                         "
                         STDC_HEADERS
                       )

CHECK_INCLUDE_FILE(sys/wait.h HAVE_SYS_WAIT_H)

CHECK_SYMBOL_EXISTS(DBL_EPSILON float.h BLT_DBL_EPSILON)
CHECK_TYPE_SIZE(int SIZEOF_INT)
CHECK_TYPE_SIZE(long SIZEOF_LONG)
CHECK_TYPE_SIZE("long long" SIZEOF_LONG_LONG)
SET(SIZEOF_VOID_P ${CMAKE_SIZEOF_VOID_P})

SET(CMAKE_REQUIRED_INCLUDES "math.h")
SET(CMAKE_REQUIRED_INCLUDES "stdlib.h")

CHECK_FUNCTION_EXISTS(drand48 HAVE_DRAND48)
IF(NOT HAVE_DRAND48)
  SET(NEED_DECL_DRAND48 1)
ENDIF(NOT HAVE_DRAND48)

CHECK_FUNCTION_EXISTS(srand48 HAVE_SRAND48)
IF(NOT HAVE_SRAND48)
  SET(NEED_DECL_SRAND48 1)
ENDIF(NOT HAVE_SRAND48)

CHECK_FUNCTION_EXISTS(strdup HAVE_STRDUP)
IF(NOT HAVE_STRDUP)
  SET(NEED_DECL_STRDUP 1)
ENDIF(NOT HAVE_STRDUP)

CHECK_FUNCTION_EXISTS(j1 HAVE_J1)
IF(NOT HAVE_J1)
  SET(NEED_DECL_J1 1)
ENDIF(NOT HAVE_J1)

CHECK_FUNCTION_EXISTS(isinfine HAVE_ISINFITE)
CHECK_FUNCTION_EXISTS(finite HAVE_INFITE)
CHECK_FUNCTION_EXISTS(isnan HAVE_ISNAN)

CHECK_FUNCTION_EXISTS(strcasecmp HAVE_STRCASECMP)
CHECK_FUNCTION_EXISTS(strncasecmp HAVE_STRNCASECMP)

CHECK_INCLUDE_FILE(ctype.h HAVE_CTYPE_H)
CHECK_INCLUDE_FILE(errno.h HAVE_ERRNO_H)
CHECK_INCLUDE_FILE(float.h HAVE_FLOAT_H)
CHECK_INCLUDE_FILE(ieeefp.h HAVE_IEEEFP_H)

CHECK_INCLUDE_FILE(inttypes.h HAVE_INTTYPES_H)
IF (NOT HAVE_INTTYPES_H)
  SET(HAVE_INTTYPES_H 0)
ENDIF (NOT HAVE_INTTYPES_H)

CHECK_INCLUDE_FILE(jpeglib.h HAVE_JPEGLIB_H)
CHECK_INCLUDE_FILE(limits.h HAVE_LIMITS_H)
CHECK_INCLUDE_FILE(math.h HAVE_MATH_H)
CHECK_INCLUDE_FILE(malloc.h HAVE_MALLOC_H)
CHECK_INCLUDE_FILE(memory.h HAVE_MEMORY_H)
CHECK_INCLUDE_FILE(setjmp.h HAVE_SETJMP_H)
CHECK_INCLUDE_FILE(stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILE(string.h HAVE_STRING_H)
CHECK_INCLUDE_FILE(unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILE(waitflags.h HAVE_WAITFLAGS_H)

IF(NOT HAVE_INTTYPES_H)
  SET(HAVE_INTTYPES_H 0)
ENDIF(NOT HAVE_INTTYPES_H)

CHECK_INCLUDE_FILE(sys/param.h HAVE_SYS_PARAM_H)
CHECK_INCLUDE_FILE(sys/time.h HAVE_SYS_TIME_H)

ADD_DEFINITIONS(-DCFG_INSTALL_LIBDIR="${CMAKE_INSTALL_PREFIX}")
ADD_DEFINITIONS(-DCFG_INSTALL_BINDIR="${CMAKE_INSTALL_PREFIX}")
ADD_DEFINITIONS(-DCFG_INSTALL_DOCDIR="${CMAKE_INSTALL_PREFIX}")

ADD_DEFINITIONS(-DCFG_RUNTIME_LIBDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
ADD_DEFINITIONS(-DCFG_RUNTIME_BINDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
ADD_DEFINITIONS(-DCFG_RUNTIME_SCRDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/library")
ADD_DEFINITIONS(-DCFG_RUNTIME_INCDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
ADD_DEFINITIONS(-DCFG_RUNTIME_DOCDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

ADD_DEFINITIONS(-DBLT_LIBRARY="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/BLT${BLT_VERSION}")
ADD_DEFINITIONS(-DBLT_PACKAGE_PATH="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
SET(HAVE_LIBM 1)

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

  # MS Windows Input Method Manager
  FIND_LIBRARY(IMM_LIBRARY imm32 PATHS ${MS_SDK_REG_PATH})
  # Common GUI components
  FIND_LIBRARY(COMCTL_LIBRARY comctl32 PATHS ${MS_SDK_REG_PATH})
  # Windows socket 2 library
  FIND_LIBRARY(WS2_LIBRARY ws2_32 PATHS ${MS_SDK_REG_PATH})
  # Windows Windows Graphical Device Interface library
  FIND_LIBRARY(GDI_LIBRARY gdi32 PATHS ${MS_SDK_REG_PATH})
ENDIF(WIN32)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/bltConfig.h.in ${CMAKE_CURRENT_BINARY_DIR}/bltConfig.h)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/bltHash.h.in ${CMAKE_CURRENT_BINARY_DIR}/bltHash.h)

#################################################################
# Set top level properties

SET(SCI_BLT_LIBRARY blt CACHE INTERNAL "BLT library name" FORCE)
SET(SCI_BLTSTUB_LIBRARY bltstub CACHE INTERNAL "BLT library name" FORCE)

SET(SCI_BLT_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/src;${CMAKE_CURRENT_BINARY_DIR}")

#################################################################
# Building library

ADD_LIBRARY(${SCI_BLT_LIBRARY} ${BLT_SRCS})
SET_TARGET_PROPERTIES(${SCI_BLT_LIBRARY}
                      PROPERTIES
                      OUTPUT_NAME BLT.${BLT_VERSION}
                     )
IF(ENABLE_PIC)
  SET_TARGET_PROPERTIES(${SCI_BLT_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fPIC")
ENDIF(ENABLE_PIC)

IF(WIN32)
  SET_TARGET_PROPERTIES(${SCI_BLT_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fp:strict -Dinline=__inline "
                          DEFINE_SYMBOL _USRDLL
                          DEFINE_SYMBOL _WINDLL
                       )

  IF(MSVC)
    SET_TARGET_PROPERTIES(${SCI_BLT_LIBRARY}
                          PROPERTIES
                            COMPILE_DEFINITIONS_RELEASE "_CRT_SECURE_NO_WARNINGS=1"
                            COMPILE_DEFINITIONS_MINSIZEREL "_CRT_SECURE_NO_WARNINGS=1"
                         )
  ENDIF(MSVC)
ENDIF(WIN32)

IF(UNIX) 
  TARGET_LINK_LIBRARIES(${SCI_BLT_LIBRARY}
                        ${SCI_TK_LIBRARY}
                        ${SCI_TCL_LIBRARY}
                        ${X11_LIBRARIES}
                        ${X11_Xft_LIB})
  IF(X11_Xscreensaver_FOUND)
    TARGET_LINK_LIBRARIES(${SCI_BLT_LIBRARY}
                            ${X11_Xscreensaver_LIB}
                         )
  ENDIF(X11_Xscreensaver_FOUND)
ELSE(UNIX) 
  IF(WIN32)
    TARGET_LINK_LIBRARIES(${SCI_BLT_LIBRARY}
                          ${SCI_TK_LIBRARY}
                          ${SCI_TCL_LIBRARY}
                          ${SCI_TCLSTUB_LIBRARY}
                          ${IMM_LIBRARY}
                          ${COMCTL_LIBRARY}
                          ${WS2_LIBRARY}
                          ${GDI_LIBRARY})
  ENDIF(WIN32)
ENDIF(UNIX) 


#################################################################
# Building directories with tcl files

FILE(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/BLT${BLT_VERSION})

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/library/tclIndex
               ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/BLT${BLT_VERSION}/tclIndex
               COPYONLY)

FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/*.tcl")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/BLT${BLT_VERSION}/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/*.xbm")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/BLT${BLT_VERSION}/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/*.cur")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/BLT${BLT_VERSION}/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/*.pro")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/BLT${BLT_VERSION}/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

SET(PACKAGE_VERSION "${BLT_VERSION}")
GET_TARGET_PROPERTY(PKG_LIB_FILE ${SCI_BLT_LIBRARY} OUTPUT_NAME)

IF(BUILD_SHARED_LIBS)
  SET(PKG_LIB_FILE "${CMAKE_SHARED_LIBRARY_PREFIX}${PKG_LIB_FILE}${CMAKE_SHARED_LIBRARY_SUFFIX}")
ELSE(BUILD_SHARED_LIBS)
  SET(PKG_LIB_FILE "${CMAKE_STATIC_LIBRARY_PREFIX}${PKG_LIB_FILE}${CMAKE_STATIC_LIBRARY_SUFFIX}")
ENDIF(BUILD_SHARED_LIBS)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/library/pkgIndex.tcl.in
               ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/BLT${BLT_VERSION}/pkgIndex.tcl)
