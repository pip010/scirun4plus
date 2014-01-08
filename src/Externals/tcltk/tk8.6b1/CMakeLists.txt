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
SET(TK_VERSION "8.6" CACHE INTERNAL "Tk version")
SET(TK_MAJOR_VERSION "8")
SET(TK_MINOR_VERSION "6")
SET(TK_PATCH_LEVEL "b1")

#################################################################
# Setting up include files and source files

SET(TK_GENERIC_SRCS
      generic/tk3d.c 
      generic/tkArgv.c
      generic/tkAtom.c 
      generic/tkBind.c
      generic/tkBitmap.c 
      generic/tkBusy.c 
      generic/tkClipboard.c
      generic/tkCmds.c 
      generic/tkColor.c
      generic/tkConfig.c 
      generic/tkCursor.c
      generic/tkError.c 
      generic/tkEvent.c
      generic/tkFocus.c
      generic/tkFont.c
      generic/tkGet.c
      generic/tkGC.c
      generic/tkGeometry.c
      generic/tkGrab.c
      generic/tkGrid.c
      generic/tkConsole.c
      generic/tkMain.c
      generic/tkOption.c
      generic/tkPack.c
      generic/tkPlace.c
      generic/tkSelect.c
      generic/tkStyle.c
      generic/tkUndo.c
      generic/tkUtil.c
      generic/tkVisual.c
      generic/tkWindow.c
      generic/tkButton.c
      generic/tkObj.c
      generic/tkEntry.c
      generic/tkFrame.c
      generic/tkListbox.c
      generic/tkMenu.c
      generic/tkMenubutton.c
      generic/tkMenuDraw.c
      generic/tkMessage.c
      generic/tkPanedWindow.c
      generic/tkScrollbar.c
      generic/tkCanvas.c
      generic/tkCanvArc.c
      generic/tkCanvBmap.c
      generic/tkCanvImg.c
      generic/tkCanvLine.c
      generic/tkCanvPoly.c
      generic/tkCanvPs.c
      generic/tkCanvText.c
      generic/tkCanvUtil.c
      generic/tkCanvWind.c
      generic/tkRectOval.c
      generic/tkScale.c
      generic/tkTrig.c
      generic/tkImage.c
      generic/tkImgBmap.c
      generic/tkImgGIF.c
      generic/tkImgPPM.c
      generic/tkImgPhoto.c
      generic/tkImgPhInstance.c
      generic/tkMacWinMenu.c
      generic/tkText.c
      generic/tkTextBTree.c
      generic/tkTextDisp.c
      generic/tkTextImage.c
      generic/tkTextIndex.c
      generic/tkTextMark.c
      generic/tkTextTag.c
      generic/tkTextWind.c
      generic/tkOldConfig.c
      generic/tkStubInit.c
   )

SET(TK_STUB_SRCS 
      generic/tkStubLib.c
   )

SET(TTK_SRCS
      generic/ttk/ttkBlink.c
      generic/ttk/ttkButton.c
      generic/ttk/ttkCache.c
      generic/ttk/ttkClamTheme.c
      generic/ttk/ttkClassicTheme.c
      generic/ttk/ttkDefaultTheme.c
      generic/ttk/ttkElements.c
      generic/ttk/ttkEntry.c
      generic/ttk/ttkFrame.c
      generic/ttk/ttkImage.c
      generic/ttk/ttkInit.c
      generic/ttk/ttkLabel.c
      generic/ttk/ttkLayout.c
      generic/ttk/ttkManager.c
      generic/ttk/ttkNotebook.c
      generic/ttk/ttkPanedwindow.c
      generic/ttk/ttkProgress.c
      generic/ttk/ttkScale.c
      generic/ttk/ttkScrollbar.c
      generic/ttk/ttkScroll.c
      generic/ttk/ttkSeparator.c
      generic/ttk/ttkState.c
      generic/ttk/ttkTagSet.c
      generic/ttk/ttkTheme.c
      generic/ttk/ttkTrace.c
      generic/ttk/ttkTrack.c
      generic/ttk/ttkTreeview.c
      generic/ttk/ttkWidget.c
      generic/ttk/ttkStubInit.c
   )

SET(TTK_STUB_SRCS
      generic/ttk/ttkStubLib.c
   )

SET(TK_UNIX_SRCS
      generic/tkSquare.c
      generic/ttk/ttkSquare.c
      unix/tkAppInit.c
      unix/tkUnix.c
      unix/tkUnix3d.c
      unix/tkUnixButton.c
      unix/tkUnixColor.c
      unix/tkUnixConfig.c
      unix/tkUnixCursor.c
      unix/tkUnixDraw.c
      unix/tkUnixEmbed.c
      unix/tkUnixEvent.c
      unix/tkUnixFocus.c
#unix/tkUnixRFont.c
      unix/tkUnixFont.c
      unix/tkUnixInit.c
      unix/tkUnixKey.c
      unix/tkUnixMenu.c
      unix/tkUnixMenubu.c
      unix/tkUnixScale.c
      unix/tkUnixScrlbr.c
      unix/tkUnixSelect.c
      unix/tkUnixSend.c 
      unix/tkUnixWm.c
      unix/tkUnixXId.c
   )

SET(TK_WIN_SRCS
      generic/tkFileFilter.c
      generic/tkImgUtil.c
      generic/tkPointer.c
      win/stubs.c
      win/tkWin32Dll.c
      win/tkWin3d.c
      win/tkWinButton.c
      win/tkWinClipboard.c
      win/tkWinColor.c
      win/tkWinConfig.c
      win/tkWinCursor.c
      win/tkWinDialog.c
      win/tkWinDraw.c
      win/tkWinEmbed.c
      win/tkWinFont.c
      win/tkWinImage.c
      win/tkWinInit.c
      win/tkWinKey.c
      win/tkWinMenu.c
      win/tkWinPixmap.c
      win/tkWinPointer.c
      win/tkWinRegion.c
      win/tkWinScrlbr.c
      win/tkWinSend.c
      win/tkWinSendCom.c
      win/tkWinWindow.c
      win/tkWinWm.c
      win/tkWinX.c
      win/ttkWinMonitor.c
      win/ttkWinTheme.c
      win/ttkWinXPTheme.c
      unix/tkUnixScale.c
      unix/tkUnixMenubu.c
   )

SET(TK_XLIB_SRCS
      xlib/xcolors.c
      xlib/xdraw.c
      xlib/xgc.c
      xlib/ximage.c
      xlib/xutil.c
   )

IF (UNIX)
  # X11 version (so we can run 64bit on OSX as well)
  SET(TK_SRCS ${TK_UNIX_SRCS} ${TTK_SRCS} ${TK_GENERIC_SRCS})
  FIND_PACKAGE( X11 )

  IF(NOT X11_FOUND)
    MESSAGE(ERROR "X11 is required to build Tk on Unix and Unix-like platforms.")
  ENDIF(NOT X11_FOUND)

  INCLUDE_DIRECTORIES(
    ${X11_INCLUDE_DIR}
    unix
    generic
    bitmaps
    macosx
    ${SCI_TCL_INCLUDE}
  )
ELSE(UNIX)
  # Windows version

  IF(WIN32)
    SET(TK_SRCS
          ${TK_WIN_SRCS}
          ${TTK_SRCS}
          ${TK_GENERIC_SRCS}
          ${TK_XLIB_SRCS}
       )

    IF(BUILD_SHARED_LIBS)
      SET(RC_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/win/rc
          CACHE INTERNAL
          "Copy windows resources to project binary directory.")
      FILE(MAKE_DIRECTORY ${RC_BINARY_DIR})
      FILE(GLOB win32_rc_files
                      "${CMAKE_CURRENT_SOURCE_DIR}/win/rc/*.bmp"
                      "${CMAKE_CURRENT_SOURCE_DIR}/win/rc/*.cur"
                      "${CMAKE_CURRENT_SOURCE_DIR}/win/rc/*.ico"
                      "${CMAKE_CURRENT_SOURCE_DIR}/win/rc/*.rc")

      FOREACH(rc_file ${win32_rc_files})
        GET_FILENAME_COMPONENT(rc_name ${rc_file} NAME)
        CONFIGURE_FILE(${rc_file}
                       ${RC_BINARY_DIR}/${rc_name}
                       COPYONLY)
      ENDFOREACH(rc_file ${win32_rc_files})

      SET(TK_SRCS ${TK_SRCS} ${RC_BINARY_DIR}/tk.rc)

      # Default machine type is IX86.
      # MACHINE=(IX86|IA64|AMD64|ALPHA)
      IF(CMAKE_SIZEOF_VOID_P MATCHES 8)
        SET(MACHINE "AMD64")
      ELSE(CMAKE_SIZEOF_VOID_P MATCHES 8)
        SET(MACHINE "IX86")
      ENDIF(CMAKE_SIZEOF_VOID_P MATCHES 8)

      SET(TK_WIN_VERSION "${TK_VERSION}.0.0")

      CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/win/wish.exe.manifest.in
                     ${RC_BINARY_DIR}/wish.exe.manifest)
    ENDIF(BUILD_SHARED_LIBS)

    INCLUDE_DIRECTORIES(
      win
      generic
      bitmaps
      xlib
      ${SCI_TCL_INCLUDE}
    )

  ENDIF(WIN32)
ENDIF(UNIX)

#################################################################
# Setting up definitions

ADD_DEFINITIONS(-DPACKAGE_NAME="tk")
ADD_DEFINITIONS(-DPACKAGE_TARNAME="tk")
ADD_DEFINITIONS(-DPACKAGE_VERSION="${TK_VERSION}")
ADD_DEFINITIONS(-DPACKAGE_STRING="tk ${TK_VERSION}")
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
ADD_DEFINITIONS(-DSTDC_HEADERS=${STDC_HEADERS})

CHECK_INCLUDE_FILE(sys/types.h HAVE_SYS_TYPES)
IF(HAVE_SYS_TYPES)
  ADD_DEFINITIONS(-DHAVE_SYS_TYPES_H=${HAVE_SYS_TYPES})
ENDIF(HAVE_SYS_TYPES)

CHECK_INCLUDE_FILE(sys/param.h HAVE_SYS_PARAM)
IF(HAVE_SYS_PARAM)
  ADD_DEFINITIONS(-DHAVE_SYS_PARAM_H=${HAVE_SYS_PARAM})
ENDIF(HAVE_SYS_PARAM)

CHECK_INCLUDE_FILE(sys/stat.h HAVE_SYS_STAT)
IF(HAVE_SYS_STAT)
  ADD_DEFINITIONS(-DHAVE_SYS_STAT_H=${HAVE_SYS_STAT})
ENDIF(HAVE_SYS_STAT)

CHECK_INCLUDE_FILE(sys/time.h HAVE_SYS_TIME)
IF(HAVE_SYS_TIME)
  ADD_DEFINITIONS(-DHAVE_SYS_TIME_H=${HAVE_SYS_TIME})
  ADD_DEFINITIONS(-DTIME_WITH_SYS_TIME_H=${HAVE_SYS_TIME})
ENDIF(HAVE_SYS_TIME)

CHECK_INCLUDE_FILE(stdlib.h HAVE_STDLIB)
IF(HAVE_STDLIB)
  ADD_DEFINITIONS(-DHAVE_STDLIB_H=${HAVE_STDLIB})
ENDIF(HAVE_STDLIB)

CHECK_INCLUDE_FILE(memory.h HAVE_MEMORY)
IF(HAVE_MEMORY)
  ADD_DEFINITIONS(-DHAVE_MEMORY_H=${HAVE_MEMORY})
ENDIF(HAVE_MEMORY)

CHECK_INCLUDE_FILE(string.h HAVE_STRING)
IF(HAVE_STRING)
  ADD_DEFINITIONS(-DHAVE_STRING_H=${HAVE_STRING})
ENDIF(HAVE_STRING)

CHECK_INCLUDE_FILE(strings.h HAVE_STRINGS)
IF(HAVE_STRINGS)
  ADD_DEFINITIONS(-DHAVE_STRINGS_H=${HAVE_STRINGS})
ENDIF(HAVE_STRINGS)

CHECK_INCLUDE_FILE(inttypes.h HAVE_INTTYPES)
IF(HAVE_INTTYPES)
  ADD_DEFINITIONS(-DHAVE_INTTYPES_H=${HAVE_INTTYPES})
ENDIF(HAVE_INTTYPES)

CHECK_INCLUDE_FILE(stdint.h HAVE_STDINT)
IF(HAVE_STDINT)
  ADD_DEFINITIONS(-DHAVE_STDINT_H=${HAVE_STDINT})
ENDIF(HAVE_STDINT)

CHECK_INCLUDE_FILE(unistd.h HAVE_UNISTD)
IF(HAVE_UNISTD)
  ADD_DEFINITIONS(-DHAVE_UNISTD_H=${HAVE_UNISTD})
ENDIF(HAVE_UNISTD)

CHECK_INCLUDE_FILE(limits.h HAVE_LIMITS)
IF(HAVE_LIMITS)
  ADD_DEFINITIONS(-DHAVE_LIMITS_H=${HAVE_LIMITS})
ENDIF(HAVE_LIMITS)

IF(UNIX)
  ADD_DEFINITIONS(-DMODULE_SCOPE=extern\ __attribute__\(\(__visibility__\(\"hidden\"\)\)\))
ENDIF(UNIX)

IF(HAVE_SYS_STAT_H)
  CHECK_C_SOURCE_COMPILES(
                           "
                             #include <sys/stat.h>

                             int main(int argc, char* argv[])
                             {
                               struct stat64 ss;
                               return 0;
                             }
                           "
                           HAVE_STRUCT_STAT64
                         )
  IF(HAVE_STRUCT_STAT64)
    ADD_DEFINITIONS(-DHAVE_STRUCT_STAT64=1)
  ENDIF(HAVE_STRUCT_STAT64)
ENDIF(HAVE_SYS_STAT_H)

ADD_DEFINITIONS(-DTCL_WIDE_INT_TYPE=long\ long)
ADD_DEFINITIONS(-DHAVE_INTPTR_T=1)
ADD_DEFINITIONS(-DHAVE_UINTPTR_T=1)
ADD_DEFINITIONS(-DTCL_CFG_OPTIMIZED=1)

CHECK_C_SOURCE_COMPILES(
                         "
                           #include <pwd.h>
                           int main(int argc, char* argv[])
                           {
                             struct passwd pwd;
                             pwd.pw_gecos;
                             return 0;
                           }
                         "
                         HAVE_PW_GECOS
                       )
IF(HAVE_PW_GECOS)
  ADD_DEFINITIONS(-DHAVE_PW_GECOS=${HAVE_PW_GECOS}) 
ENDIF(HAVE_PW_GECOS)

IF(X11_Xft_FOUND)
#  ADD_DEFINITIONS(-DHAVE_LIBXFT=1)
#  ADD_DEFINITIONS(-DHAVE_XFT=1)
ENDIF(X11_Xft_FOUND)

IF(X11_Xscreensaver_FOUND)
  ADD_DEFINITIONS(-DHAVE_XSS=1)
ENDIF(X11_Xscreensaver_FOUND)

IF(APPLE)
  ADD_DEFINITIONS(-DMAC_OSX_TCL=1)
  ADD_DEFINITIONS(-DHAVE_COREFOUNDATION=1)
ENDIF(APPLE)

ADD_DEFINITIONS(-DCFG_INSTALL_LIBDIR="${CMAKE_INSTALL_PREFIX}")
ADD_DEFINITIONS(-DCFG_INSTALL_BINDIR="${CMAKE_INSTALL_PREFIX}")
ADD_DEFINITIONS(-DCFG_INSTALL_DOCDIR="${CMAKE_INSTALL_PREFIX}")

ADD_DEFINITIONS(-DCFG_RUNTIME_LIBDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
ADD_DEFINITIONS(-DCFG_RUNTIME_BINDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
ADD_DEFINITIONS(-DCFG_RUNTIME_SCRDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/library")
ADD_DEFINITIONS(-DCFG_RUNTIME_INCDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
ADD_DEFINITIONS(-DCFG_RUNTIME_DOCDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

ADD_DEFINITIONS(-DTK_LIBRARY="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tk${TK_VERSION}")
ADD_DEFINITIONS(-DTK_PACKAGE_PATH="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

IF(WIN32)
  # MS Windows Input Method Manager
  FIND_LIBRARY(IMM_LIBRARY imm32 PATHS ${MS_SDK_REG_PATH})
  # Common GUI components
  FIND_LIBRARY(COMCTL_LIBRARY comctl32 PATHS ${MS_SDK_REG_PATH})
  # Windows socket 2 library
  FIND_LIBRARY(WS2_LIBRARY ws2_32 PATHS ${MS_SDK_REG_PATH})
  # Windows Windows Graphical Device Interface library
  FIND_LIBRARY(GDI_LIBRARY gdi32 PATHS ${MS_SDK_REG_PATH})
ENDIF(WIN32)

#################################################################
# Set top level properties

GET_DIRECTORY_PROPERTY(TK_DEFINITIONS COMPILE_DEFINITIONS)

SET(SCI_TK_LIBRARY tk CACHE INTERNAL "TK library name" FORCE)
SET(SCI_TKSTUB_LIBRARY tkstub CACHE INTERNAL "TK library name" FORCE)

IF(UNIX)
  SET(SCI_TK_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/generic;${CMAKE_CURRENT_SOURCE_DIR}/unix;${X11_INCLUDE_DIR}"
      CACHE INTERNAL "TK include directories"
      FORCE)
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
      STRING(REGEX REPLACE "-O3" "-O0" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
    ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
ELSE(UNIX)
  IF(WIN32)
    SET(SCI_TK_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/generic;${CMAKE_CURRENT_SOURCE_DIR}/unix;${CMAKE_CURRENT_SOURCE_DIR}/win;${CMAKE_CURRENT_SOURCE_DIR}/xlib"
        CACHE INTERNAL "TK include directories"
        FORCE)

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
  ENDIF(WIN32)
ENDIF(UNIX)

SET(SCI_TK_DEFINITIONS "${TK_DEFINITIONS}" CACHE INTERNAL "Definitions used for tk" FORCE)

#################################################################
# Building library

ADD_LIBRARY(${SCI_TK_LIBRARY} ${TK_SRCS})
SET_TARGET_PROPERTIES(${SCI_TK_LIBRARY}
                      PROPERTIES
                        OUTPUT_NAME tk.${TK_VERSION}
                        COMPILE_DEFINITIONS_DEBUG "TCL_CFG_DEBUG=1"
                     )

IF(WIN32)
  SET_TARGET_PROPERTIES(${SCI_TK_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fp:strict -Dinline=__inline"
                          COMPILE_DEFINITIONS "BUILD_tk;BUILD_ttk;SUPPORT_CONFIG_EMBEDDED;_MBCS"
                          DEFINE_SYMBOL _USRDLL
                          DEFINE_SYMBOL _WINDLL
                     )
  IF(MSVC)
    SET_TARGET_PROPERTIES(${SCI_TK_LIBRARY}
                          PROPERTIES
                            COMPILE_DEFINITIONS_RELEASE "_CRT_SECURE_NO_WARNINGS=1"
                            COMPILE_DEFINITIONS_MINSIZEREL "_CRT_SECURE_NO_WARNINGS=1"
                         )
  ENDIF(MSVC)
ENDIF(WIN32)

IF(UNIX) 
  TARGET_LINK_LIBRARIES(${SCI_TK_LIBRARY}
                        ${SCI_TCL_LIBRARY}
                        ${X11_LIBRARIES}
                        ${X11_Xft_LIB})
  IF(X11_Xscreensaver_FOUND)
    TARGET_LINK_LIBRARIES(${SCI_TK_LIBRARY}
                            ${X11_Xscreensaver_LIB}
                         )
  ENDIF(X11_Xscreensaver_FOUND)
ELSE(UNIX) 
  IF(WIN32)
    TARGET_LINK_LIBRARIES(${SCI_TK_LIBRARY}
                          ${SCI_TCL_LIBRARY}
                          ${SCI_TCLSTUB_LIBRARY}
                          ${IMM_LIBRARY}
                          ${COMCTL_LIBRARY}
                          ${WS2_LIBRARY}
                          ${GDI_LIBRARY})
  ENDIF(WIN32)
ENDIF(UNIX) 

IF(APPLE)
  SET_TARGET_PROPERTIES(${SCI_TK_LIBRARY} PROPERTIES LINK_FLAGS  "-framework CoreFoundation")
ENDIF(APPLE)

ADD_LIBRARY(${SCI_TKSTUB_LIBRARY} STATIC ${TK_STUB_SRCS})
SET_TARGET_PROPERTIES(${SCI_TKSTUB_LIBRARY}
                      PROPERTIES
                        OUTPUT_NAME tkstub.${TK_VERSION}
                        COMPILE_DEFINITIONS "STATIC_BUILD"
                        COMPILE_DEFINITIONS_DEBUG "TCL_CFG_DEBUG=1"
                     )
IF(UNIX AND NOT APPLE AND NOT CYGWIN)
  SET_TARGET_PROPERTIES(${SCI_TKSTUB_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fPIC"
                       )
ENDIF(UNIX AND NOT APPLE AND NOT CYGWIN)

IF(ENABLE_PIC)
  SET_TARGET_PROPERTIES(${SCI_TK_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fPIC")
ENDIF(ENABLE_PIC)

#################################################################
# Building directories with tcl files

FILE(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tk${TK_VERSION})

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/library/tclIndex
               ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tk${TK_VERSION}/tclIndex
               COPYONLY)

FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/*.tcl")

FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tk${TK_VERSION}/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

FILE(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tk${TK_VERSION}/ttk)
FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/ttk/*.tcl")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tk${TK_VERSION}/ttk/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

FILE(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tk${TK_VERSION}/msgs)
FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/msgs/*.msg")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tk${TK_VERSION}/msgs/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

FILE(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tk${TK_VERSION}/images)
FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/images/*.gif")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tk${TK_VERSION}/images/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/images/*.eps")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tk${TK_VERSION}/images/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})


SET(PACKAGE_VERSION "${TK_VERSION}")
GET_TARGET_PROPERTY(PKG_LIB_FILE ${SCI_TK_LIBRARY} OUTPUT_NAME)

IF(BUILD_SHARED_LIBS)
  SET(PKG_LIB_FILE "${CMAKE_SHARED_LIBRARY_PREFIX}${PKG_LIB_FILE}${CMAKE_SHARED_LIBRARY_SUFFIX}")
ELSE(BUILD_SHARED_LIBS)
  SET(PKG_LIB_FILE "${CMAKE_STATIC_LIBRARY_PREFIX}${PKG_LIB_FILE}${CMAKE_STATIC_LIBRARY_SUFFIX}")
ENDIF(BUILD_SHARED_LIBS)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/library/pkgIndex.tcl.in
               ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tk${TK_VERSION}/pkgIndex.tcl)
