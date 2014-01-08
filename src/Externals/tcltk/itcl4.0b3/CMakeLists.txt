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
SET(ITCL_VERSION "4.0b3" CACHE INTERNAL "ITCL version")
SET(ITCL_MAJOR_VERSION "4")
SET(ITCL_MINOR_VERSION "0")
SET(ITCL_PATCH_LEVEL "b3")


SET(ITCL_GENERIC_SRCS
      generic/itcl2TclOO.c 
      generic/itclBase.c
      generic/itclBuiltin.c
      generic/itclClass.c
      generic/itclCmd.c
      generic/itclEnsemble.c
      generic/itclHelpers.c
      generic/itclInfo.c
      generic/itclLinkage.c
      generic/itclMethod.c
      generic/itclObject.c
      generic/itclParse.c
      generic/itclStubs.c 
      generic/itclStubInit.c 
      generic/itclResolve.c
      generic/itclTclIntStubsFcn.c 
      generic/itclUtil.c 
      generic/itclMigrate2TclCore.c 
      generic/itclNeededFromTclOO.c
)

SET(ITCL_STUB_SRCS
      generic/itclStubLib.c
)

SET(ITCL_WIN32_SRCS
      win/dllEntryPoint.c
   )


IF(UNIX)
  # LINUX/APPLE BUILD
  SET(ITCL_SRCS ${ITCL_GENERIC_SRCS} ) 
ELSE(UNIX)
  # WINDOWS BUILD
  IF(WIN32)
    SET(ITCL_SRCS
        ${ITCL_GENERIC_SRCS}
        ${ITCL_WIN32_SRCS}
       ) 
    IF(BUILD_SHARED_LIBS)
      SET(ITCL_SRCS ${ITCL_SRCS} win/itcl.rc)
    ENDIF(BUILD_SHARED_LIBS)
  ENDIF(WIN32)
ENDIF(UNIX)

INCLUDE_DIRECTORIES(
  ${SCI_TCL_INCLUDE}
  ${SCI_TK_INCLUDE}
  generic
  win
)

#################################################################
# Setting up definitions

ADD_DEFINITIONS(-DPACKAGE_NAME="itcl")
ADD_DEFINITIONS(-DPACKAGE_TARNAME="itcl")
ADD_DEFINITIONS(-DPACKAGE_VERSION="${ITCL_VERSION}")
ADD_DEFINITIONS(-DPACKAGE_STRING="itcl ${ITCL_VERSION}")
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

CHECK_INCLUDE_FILE(sys/stat.h HAVE_SYS_STAT)
IF(HAVE_SYS_STAT)
  ADD_DEFINITIONS(-DHAVE_SYS_STAT_H=${HAVE_SYS_STAT})
ENDIF(HAVE_SYS_STAT)

CHECK_INCLUDE_FILE(sys/param.h HAVE_SYS_PARAM)
IF(HAVE_SYS_PARAM)
  ADD_DEFINITIONS(-DHAVE_SYS_PARAM_H=${HAVE_SYS_PARAM})
ENDIF(HAVE_SYS_PARAM)

CHECK_INCLUDE_FILE(stdlib.h HAVE_STDLIB)
IF(HAVE_STDLIB)
  ADD_DEFINITIONS(-DHAVE_STDLIB_H=${HAVE_STDLIB})
ENDIF(HAVE_STDLIB)

CHECK_INCLUDE_FILE(string.h HAVE_STRING)
IF(HAVE_STRING)
  ADD_DEFINITIONS(-DHAVE_STRING_H=${HAVE_STRING})
ENDIF(HAVE_STRING)

CHECK_INCLUDE_FILE(strings.h HAVE_STRINGS)
IF(HAVE_STRINGS)
  ADD_DEFINITIONS(-DHAVE_STRINGS_H=${HAVE_STRINGS})
ENDIF(HAVE_STRINGS)

CHECK_INCLUDE_FILE(memory.h HAVE_MEMORY)
IF(HAVE_MEMORY)
  ADD_DEFINITIONS(-DHAVE_MEMORY_H=${HAVE_MEMORY})
ENDIF(HAVE_MEMORY)

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

CHECK_INCLUDE_FILE(values.h HAVE_VALUES)
SET(NO_VALUES "")
IF(NOT HAVE_VALUES)
  SET(NO_VALUES 1)
ENDIF(NOT HAVE_VALUES)
ADD_DEFINITIONS(-DNO_VALUES_H=${NO_VALUES})

CHECK_INCLUDE_FILE(limits.h HAVE_LIMITS)
IF(HAVE_LIMITS)
  ADD_DEFINITIONS(-DHAVE_LIMITS_H=${HAVE_LIMITS})
ENDIF(HAVE_LIMITS)

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

IF(UNIX)
  ADD_DEFINITIONS(-DMODULE_SCOPE=extern\ __attribute__\(\(__visibility__\(\"hidden\"\)\)\))
ENDIF(UNIX)

ADD_DEFINITIONS(-DTCL_WIDE_INT_TYPE=long\ long)
ADD_DEFINITIONS(-DUSE_TCL_STUBS=1)
ADD_DEFINITIONS(-DUSE_TCLOO_STUBS=1)

ADD_DEFINITIONS(-DITCL_LIBRARY="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/itcl${ITCL_VERSION}")
ADD_DEFINITIONS(-DITCL_PACKAGE_PATH="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")


#################################################################
# Set top level properties

GET_DIRECTORY_PROPERTY(ITCL_DEFINITIONS COMPILE_DEFINITIONS)

SET(SCI_ITCL_LIBRARY itcl CACHE INTERNAL "ITCL library name" FORCE)
SET(SCI_ITCLSTUB_LIBRARY itclstub CACHE INTERNAL "ITCL library name" FORCE)
SET(SCI_ITCL_DEFINITIONS "${ITCL_DEFINITIONS}" CACHE INTERNAL "Definitions used for itcl" FORCE)
SET(SCI_ITCL_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/generic"
    CACHE INTERNAL "ITCL include directories"
    FORCE)

IF(WIN32 AND MSVC)
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
ENDIF(WIN32 AND MSVC)

#################################################################
# Building library

ADD_LIBRARY(${SCI_ITCL_LIBRARY} ${ITCL_SRCS})
TARGET_LINK_LIBRARIES(${SCI_ITCL_LIBRARY}
                      ${SCI_TCLSTUB_LIBRARY}
                      ${SCI_TKSTUB_LIBRARY})
SET_TARGET_PROPERTIES(${SCI_ITCL_LIBRARY}
                      PROPERTIES
                        OUTPUT_NAME itcl.${ITCL_VERSION})
IF(ENABLE_PIC)
  SET_TARGET_PROPERTIES(${SCI_ITCL_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fPIC")
ENDIF(ENABLE_PIC)

IF(APPLE)
  SET_TARGET_PROPERTIES(${SCI_ITCL_LIBRARY}
                        PROPERTIES
                          LINK_FLAGS "-framework CoreFoundation")
ENDIF(APPLE)

IF(WIN32)
  SET_TARGET_PROPERTIES(${SCI_ITCL_LIBRARY}
                        PROPERTIES
                          DEFINE_SYMBOL _USRDLL
                          DEFINE_SYMBOL _WINDLL
                          COMPILE_DEFINITIONS "_MBCS;TCL_THREADS=1;BUILD_itcl"
                       )
  IF(MSVC)
    SET_TARGET_PROPERTIES(${SCI_ITCL_LIBRARY}
                          PROPERTIES
                            COMPILE_DEFINITIONS_RELEASE "_CRT_SECURE_NO_WARNINGS=1"
                            COMPILE_DEFINITIONS_MINSIZEREL "_CRT_SECURE_NO_WARNINGS=1"
                         )
  ENDIF(MSVC)
ENDIF(WIN32)

ADD_LIBRARY(${SCI_ITCLSTUB_LIBRARY} STATIC ${ITCL_STUB_SRCS})
SET_TARGET_PROPERTIES(${SCI_ITCLSTUB_LIBRARY}
                      PROPERTIES
                        OUTPUT_NAME itclstub.${ITCL_VERSION}
                        COMPILE_DEFINITIONS "STATIC_BUILD"
                      )
IF(UNIX AND NOT APPLE AND NOT CYGWIN)
  SET_TARGET_PROPERTIES(${SCI_ITCLSTUB_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fPIC"
                       )
ENDIF(UNIX AND NOT APPLE AND NOT CYGWIN)

#################################################################
# Building directories with tcl files

FILE(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/itcl${ITCL_VERSION})

FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/*.tcl")

FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/itcl${ITCL_VERSION}/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

SET(PACKAGE_VERSION "${ITCL_VERSION}")
GET_TARGET_PROPERTY(PKG_LIB_FILE ${SCI_ITCL_LIBRARY} OUTPUT_NAME)

IF(BUILD_SHARED_LIBS)
  SET(PKG_LIB_FILE "${CMAKE_SHARED_LIBRARY_PREFIX}${PKG_LIB_FILE}${CMAKE_SHARED_LIBRARY_SUFFIX}")
ELSE(BUILD_SHARED_LIBS)
  SET(PKG_LIB_FILE "${CMAKE_STATIC_LIBRARY_PREFIX}${PKG_LIB_FILE}${CMAKE_STATIC_LIBRARY_SUFFIX}")
ENDIF(BUILD_SHARED_LIBS)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/pkgIndex.tcl.in
               ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/itcl${ITCL_VERSION}/pkgIndex.tcl)

