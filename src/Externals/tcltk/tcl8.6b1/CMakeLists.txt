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
SET(TCL_VERSION "8.6" CACHE INTERNAL "Tcl version")
SET(TCL_MAJOR_VERSION "8")
SET(TCL_MINOR_VERSION "6")
SET(TCL_PATCH_LEVEL "b1")

IF(UNIX)
  SET(TCLSH_SRCS
        unix/tclAppInit.c
     )
ELSE(UNIX)
  SET(TCLSH_SRCS
        win/tclAppInit.c
     )

  IF(BUILD_SHARED_LIBS)
    SET(TCLSH_SRCS ${TCLSH_SRCS} win/tclsh.rc)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(UNIX)

SET(TCL_GENERIC_SRCS
      generic/regcomp.c
      generic/regexec.c
      generic/regfree.c
      generic/regerror.c
      generic/tclAlloc.c
      generic/tclAsync.c
      generic/tclBasic.c
      generic/tclBinary.c
      generic/tclCkalloc.c
      generic/tclClock.c
      generic/tclCmdAH.c
      generic/tclCmdIL.c
      generic/tclCmdMZ.c
      generic/tclCompCmds.c
      generic/tclCompExpr.c
      generic/tclCompile.c
      generic/tclConfig.c
      generic/tclDate.c
      generic/tclDictObj.c
      generic/tclEncoding.c
      generic/tclEnv.c
      generic/tclEvent.c
      generic/tclExecute.c
      generic/tclFCmd.c
      generic/tclFileName.c
      generic/tclGet.c
      generic/tclHash.c
      generic/tclHistory.c
      generic/tclIndexObj.c
      generic/tclInterp.c
      generic/tclIO.c
      generic/tclIOCmd.c
      generic/tclIOGT.c
      generic/tclIOSock.c
      generic/tclIOUtil.c
      generic/tclIORChan.c
      generic/tclIORTrans.c
      generic/tclLink.c
      generic/tclListObj.c
      generic/tclLiteral.c
      generic/tclLoad.c
      generic/tclMain.c
      generic/tclNamesp.c
      generic/tclNotify.c
      generic/tclObj.c
      generic/tclPanic.c
      generic/tclParse.c
      generic/tclPathObj.c
      generic/tclPipe.c
      generic/tclPkg.c
      generic/tclPkgConfig.c
      generic/tclPosixStr.c
      generic/tclPreserve.c
      generic/tclProc.c
      generic/tclRegexp.c
      generic/tclResolve.c
      generic/tclResult.c
      generic/tclScan.c
      generic/tclStringObj.c
      generic/tclStrToD.c
      generic/tclStubInit.c
      generic/tclTestProcBodyObj.c
      generic/tclThread.c
      generic/tclThreadAlloc.c
      generic/tclThreadJoin.c
      generic/tclThreadStorage.c
      generic/tclTimer.c
      generic/tclTomMathInterface.c
      generic/tclTrace.c
      generic/tclUtil.c
      generic/tclUtf.c
      generic/tclUniData.c
      generic/tclVar.c
      generic/tclZlib.c
)

SET(TCL_OO_SRCS
      generic/tclOO.c
      generic/tclOOBasic.c
      generic/tclOOCall.c
      generic/tclOODefineCmds.c
      generic/tclOOInfo.c
      generic/tclOOMethod.c
      generic/tclOOStubInit.c
)

SET(TCL_STUB_SRCS
      generic/tclStubLib.c
      generic/tclOOStubLib.c
   )


SET(TCL_TOMMATH_SRCS
      libtommath/bncore.c
      libtommath/bn_reverse.c
      libtommath/bn_fast_s_mp_mul_digs.c
      libtommath/bn_fast_s_mp_sqr.c
      libtommath/bn_mp_add.c
      libtommath/bn_mp_add_d.c
      libtommath/bn_mp_and.c
      libtommath/bn_mp_clamp.c
      libtommath/bn_mp_clear.c
      libtommath/bn_mp_clear_multi.c
      libtommath/bn_mp_cmp.c
      libtommath/bn_mp_cmp_d.c
      libtommath/bn_mp_cmp_mag.c
      libtommath/bn_mp_copy.c
      libtommath/bn_mp_count_bits.c
      libtommath/bn_mp_div.c
      libtommath/bn_mp_div_d.c
      libtommath/bn_mp_div_2.c
      libtommath/bn_mp_div_2d.c
      libtommath/bn_mp_div_3.c
      libtommath/bn_mp_exch.c
      libtommath/bn_mp_expt_d.c
      libtommath/bn_mp_grow.c
      libtommath/bn_mp_init.c
      libtommath/bn_mp_init_copy.c
      libtommath/bn_mp_init_multi.c
      libtommath/bn_mp_init_set.c
      libtommath/bn_mp_init_size.c
      libtommath/bn_mp_karatsuba_mul.c
      libtommath/bn_mp_karatsuba_sqr.c
      libtommath/bn_mp_lshd.c
      libtommath/bn_mp_mod.c
      libtommath/bn_mp_mod_2d.c
      libtommath/bn_mp_mul.c
      libtommath/bn_mp_mul_2.c
      libtommath/bn_mp_mul_2d.c
      libtommath/bn_mp_mul_d.c
      libtommath/bn_mp_neg.c
      libtommath/bn_mp_or.c
      libtommath/bn_mp_radix_size.c
      libtommath/bn_mp_radix_smap.c
      libtommath/bn_mp_read_radix.c
      libtommath/bn_mp_rshd.c
      libtommath/bn_mp_set.c
      libtommath/bn_mp_shrink.c
      libtommath/bn_mp_sqr.c
      libtommath/bn_mp_sqrt.c
      libtommath/bn_mp_sub.c
      libtommath/bn_mp_sub_d.c
      libtommath/bn_mp_to_unsigned_bin.c
      libtommath/bn_mp_to_unsigned_bin_n.c
      libtommath/bn_mp_toom_mul.c
      libtommath/bn_mp_toom_sqr.c
      libtommath/bn_mp_toradix_n.c
      libtommath/bn_mp_unsigned_bin_size.c
      libtommath/bn_mp_xor.c
      libtommath/bn_mp_zero.c
      libtommath/bn_s_mp_add.c
      libtommath/bn_s_mp_mul_digs.c
      libtommath/bn_s_mp_sqr.c
      libtommath/bn_s_mp_sub.c
    )

SET(TCL_UNIX_SRCS
      unix/tclAppInit.c
      unix/tclUnixChan.c
      unix/tclUnixEvent.c
      unix/tclUnixFCmd.c
      unix/tclUnixFile.c
      unix/tclUnixPipe.c
      unix/tclUnixSock.c
      unix/tclUnixTest.c
      unix/tclUnixThrd.c
      unix/tclUnixTime.c
      unix/tclUnixInit.c
      unix/tclUnixCompat.c
    )

SET(TCL_WIN32_SRCS
      win/tclWin32Dll.c
      win/tclWinChan.c
      win/tclWinConsole.c
      win/tclWinError.c
      win/tclWinFCmd.c
      win/tclWinFile.c
      win/tclWinInit.c
      win/tclWinLoad.c
      win/tclWinPipe.c
      win/tclWinSerial.c
      win/tclWinSock.c
      win/tclWinThrd.c
      win/tclWinTime.c
   )

SET(TCL_WIN32_PIPE_SRCS
      win/stub16.c
   )

SET(TCL_WIN32_COMPAT_SRCS
      compat/fixstrtod.c
      compat/gettod.c
      compat/memcmp.c
      compat/opendir.c
      compat/strncasecmp.c
      compat/strstr.c
      compat/strtod.c
      compat/strtol.c
      compat/strtoul.c
      compat/waitpid.c
   )

IF(BUILD_SHARED_LIBS)
  SET(TCL_WIN32_SRCS ${TCL_WIN32_SRCS} win/tcl.rc)
ENDIF(BUILD_SHARED_LIBS)


IF(UNIX)
  SET(TCL_NOTIFY_SRCS
        unix/tclUnixNotfy.c
     )
ELSE(UNIX)
  IF(WIN32)
    SET(TCL_NOTIFY_SRCS
          win/tclWinNotify.c
       )
  ENDIF(WIN32)
ENDIF(UNIX)

SET(TCL_MAC_OSX_SRCS 
      macosx/tclMacOSXBundle.c
      macosx/tclMacOSXFCmd.c
      macosx/tclMacOSXNotify.c
   )

IF(UNIX)
  IF(APPLE)
    # APPLE BUILD
    SET(TCL_SRCS
          ${TCL_GENERIC_SRCS}
          ${TCL_TOMMATH_SRCS}
          ${TCL_UNIX_SRCS}
          ${TCL_OO_SRCS}
          ${TCL_STUB_SRCS}
          ${TCL_NOTIFY_SRCS}
          ${TCL_MAC_OSX_SRCS}
          unix/tclLoadDyld.c
       )
  ELSE(APPLE)
    # LINUX BUILD
    SET(TCL_SRCS
          ${TCL_GENERIC_SRCS}
          ${TCL_TOMMATH_SRCS}
          ${TCL_UNIX_SRCS}
          ${TCL_OO_SRCS}
          ${TCL_STUB_SRCS}
          ${TCL_NOTIFY_SRCS}
          unix/tclLoadDl.c
       )
  ENDIF(APPLE)
ELSE(UNIX)
    # WINDOWS BUILD
    IF(WIN32)
      SET(TCL_SRCS
            ${TCL_GENERIC_SRCS}
            ${TCL_OO_SRCS}
            ${TCL_WIN32_SRCS}
            ${TCL_NOTIFY_SRCS}
            ${TCL_TOMMATH_SRCS}
         )
    ENDIF(WIN32)
ENDIF(UNIX)

#################################################################
# Setting up include files


IF(UNIX)
  IF(APPLE)
    INCLUDE_DIRECTORIES(
                         regex
                         generic
                         libtommath
                         macosx
                         unix
                       )
  ELSE(APPLE)
    # LINUX BUILD
    INCLUDE_DIRECTORIES(
                         regex
                         generic
                         libtommath
                         unix
                         ${SCI_ZLIB_INCLUDE}
                       )
  ENDIF(APPLE)
ELSE(UNIX)
 IF(WIN32)
    INCLUDE_DIRECTORIES(
                         regex
                         generic
                         libtommath
                         win
                         ${SCI_ZLIB_INCLUDE}
                     )
  ENDIF(WIN32)
ENDIF(UNIX)


#################################################################
# Setting up definitions

ADD_DEFINITIONS(-DPACKAGE_NAME="tcl")
ADD_DEFINITIONS(-DPACKAGE_TARNAME="tcl")
ADD_DEFINITIONS(-DPACKAGE_VERSION="${TCL_VERSION}")
ADD_DEFINITIONS(-DPACKAGE_STRING="tcl ${TCL_VERSION}")
ADD_DEFINITIONS(-DPACKAGE_BUGREPORT="")

IF(UNIX)
  ADD_DEFINITIONS(-DUSE_TCLALLOC=0)
ENDIF(UNIX)

# Check for Std headers
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

CHECK_INCLUDE_FILE(sys/ioctl.h HAVE_SYS_IOCTL)
IF(HAVE_SYS_IOCTL)
  ADD_DEFINITIONS(-DHAVE_SYS_IOCTL_H=${HAVE_SYS_IOCTL})
ENDIF(HAVE_SYS_IOCTL)

CHECK_INCLUDE_FILE(time.h HAVE_TIME_H)
IF(HAVE_TIME_H)
  ADD_DEFINITIONS(-DHAVE_TIME_H=${HAVE_TIME_H})
ENDIF(HAVE_TIME_H)

CHECK_INCLUDE_FILE(sys/stat.h HAVE_SYS_STAT)
IF(HAVE_SYS_STAT)
  ADD_DEFINITIONS(-DHAVE_SYS_STAT_H=${HAVE_SYS_STAT})
ENDIF(HAVE_SYS_STAT)

CHECK_INCLUDE_FILE(sys/time.h HAVE_SYS_TIME)
IF(HAVE_SYS_TIME)
  ADD_DEFINITIONS(-DHAVE_SYS_TIME_H=${HAVE_SYS_TIME})
  ADD_DEFINITIONS(-DTIME_WITH_SYS_TIME=${HAVE_SYS_TIME})
ENDIF(HAVE_SYS_TIME)

CHECK_INCLUDE_FILE(sys/filio.h HAVE_SYS_FILIO)
IF(HAVE_SYS_FILIO)
  ADD_DEFINITIONS(-DHAVE_SYS_FILIO_H=${HAVE_SYS_FILIO})
ENDIF(HAVE_SYS_FILIO)

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

CHECK_INCLUDE_FILE(values.h HAVE_VALUES)
SET(NO_VALUES "")
IF(NOT HAVE_VALUES)
  SET(NO_VALUES 1)
ENDIF(NOT HAVE_VALUES)
ADD_DEFINITIONS(-DNO_VALUES_H=${NO_VALUES})

CHECK_INCLUDE_FILE(fts.h HAVE_FTS)
IF(HAVE_FTS)
  ADD_DEFINITIONS(-DHAVE_FTS=${HAVE_FTS})
ENDIF(HAVE_FTS)

CHECK_INCLUDE_FILE(termios.h HAVE_TERMIOS)
IF(HAVE_TERMIOS)
  ADD_DEFINITIONS(-DUSE_TERMIOS=${HAVE_TERMIOS})
ENDIF(HAVE_TERMIOS)

CHECK_INCLUDE_FILE(AvailabilityMacros.h HAVE_AVAILABITITYMAROS)
IF(HAVE_AVAILABITITYMAROS)
  ADD_DEFINITIONS(-DHAVE_AVAILABILITYMACROS_H=${HAVE_AVAILABITITYMAROS})
ENDIF(HAVE_AVAILABITITYMAROS)

SET(CMAKE_REQUIRED_INCLUDES "limits.h")
CHECK_FUNCTION_EXISTS(waitpid HAVE_WAITPID)
IF(HAVE_WAITPID)
  ADD_DEFINITIONS(-DHAVE_WAITPID=${HAVE_WAITPID})
ENDIF(HAVE_WAITPID)

CHECK_FUNCTION_EXISTS(opendir HAVE_OPENDIR)
IF(HAVE_OPENDIR)
  ADD_DEFINITIONS(-DHAVE_OPENDIR=${HAVE_OPENDIR})
ENDIF(HAVE_OPENDIR)

CHECK_FUNCTION_EXISTS(strtol HAVE_STRTOL)
IF(HAVE_STRTOL)
  ADD_DEFINITIONS(-DHAVE_STRTOL=${HAVE_STRTOL})
ENDIF(HAVE_STRTOL)

CHECK_FUNCTION_EXISTS(getcwd HAVE_GETCWD)
IF(HAVE_GETCWD)
  ADD_DEFINITIONS(-DHAVE_GETCWD=${HAVE_GETCWD})
ENDIF(HAVE_GETCWD)

CHECK_FUNCTION_EXISTS(getaddrinfo HAVE_GET_ADDRINFO)
IF(HAVE_GETADDRINFO)
  ADD_DEFINITIONS(-DHAVE_GETADDRINFO=${HAVE_GETADDRINFO})
ENDIF(HAVE_GETADDRINFO)

ADD_DEFINITIONS(-DTCL_CFG_OPTIMIZED=1)
ADD_DEFINITIONS(-DTCL_TOMMATH=1)
ADD_DEFINITIONS(-DMP_PREC=4)
ADD_DEFINITIONS(-DTCL_WIDE_INT_TYPE=long\ long)
ADD_DEFINITIONS(-DTCL_CFGVAL_ENCODING="iso8859-1")
ADD_DEFINITIONS(-DHAVE_STRUCT_STAT_ST_BLKSIZE=1)
ADD_DEFINITIONS(-DHAVE_ST_BLKSIZE=1)
ADD_DEFINITIONS(-DHAVE_INTPTR_T=1)
ADD_DEFINITIONS(-DHAVE_UINTPTR_T=1) 
ADD_DEFINITIONS(-DHAVE_SIGNED_CHAR=1)

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

IF(UNIX)
  ADD_DEFINITIONS(-DHAVE_PTHREAD_ATFORK=1)
ENDIF(UNIX)

# We know that on the newer OSX platforms we need to include these
# Hence we do not explicitly test for them
IF(APPLE)
  ADD_DEFINITIONS(-DMAC_OSX_TCL=1 
                  -DHAVE_COREFOUNDATION=1 
                  -DHAVE_GETATTRLIST=1 
                  -DHAVE_COPYFILE_H=1 
                  -DHAVE_COPYFILE=1 
                  -DUSE_VFORK=1
                  -DTCL_WIDE_CLICKS=1
                  -DTCL_DEFAULT_ENCODING="utf-8"
                  -DTCL_LOAD_FROM_MEMORY=1
                  -DHAVE_WEAK_IMPORT=1
                  )
ENDIF(APPLE)



ADD_DEFINITIONS(-DCFG_INSTALL_LIBDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
ADD_DEFINITIONS(-DCFG_INSTALL_BINDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
ADD_DEFINITIONS(-DCFG_INSTALL_SCRDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}")
ADD_DEFINITIONS(-DCFG_INSTALL_INCDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
ADD_DEFINITIONS(-DCFG_INSTALL_DOCDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

ADD_DEFINITIONS(-DCFG_RUNTIME_LIBDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
ADD_DEFINITIONS(-DCFG_RUNTIME_BINDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
ADD_DEFINITIONS(-DCFG_RUNTIME_SCRDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}")
ADD_DEFINITIONS(-DCFG_RUNTIME_INCDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
ADD_DEFINITIONS(-DCFG_RUNTIME_DOCDIR="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

ADD_DEFINITIONS(-DTCL_LIBRARY="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}")
ADD_DEFINITIONS(-DTCL_PACKAGE_PATH="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

#################################################################
# Building directories with tcl files

FILE(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION})

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/library/tclIndex
               ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}/tclIndex
               COPYONLY)
FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/*.tcl")

FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

FILE(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}/encoding)
FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/encoding/*.enc")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}/encoding/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

FILE(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}/http1.0)
FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/http1.0/*.tcl")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}/http1.0/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

FILE(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}/opt0.4)
FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/opt0.4/*.tcl")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}/opt0.4/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})


FILE(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}/msgs)
FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/msgs/*.msg")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}/msgs/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})

FILE(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}/msgcat)
FILE(GLOB tcl_files "${CMAKE_CURRENT_SOURCE_DIR}/library/msgcat/*.tcl")
FOREACH(tcl_file ${tcl_files})
  GET_FILENAME_COMPONENT(tcl_name ${tcl_file} NAME)
  CONFIGURE_FILE(${tcl_file}
                 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tcl${TCL_VERSION}/msgcat/${tcl_name}
                 COPYONLY)
ENDFOREACH(tcl_file ${tcl_files})


#################################################################
# Set top level properties

GET_DIRECTORY_PROPERTY(TCL_DEFINITIONS COMPILE_DEFINITIONS)

SET(SCI_TCL_LIBRARY tcl CACHE INTERNAL "TCL library name" FORCE)
SET(SCI_TCLSTUB_LIBRARY tclstub CACHE INTERNAL "TCL library name" FORCE)
SET(SCI_TCL_DEFINITIONS "${TCL_DEFINITIONS}" CACHE INTERNAL "Definitions used for tcl" FORCE)

IF(UNIX)
  SET(SCI_TCL_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/generic;${CMAKE_CURRENT_SOURCE_DIR}/unix"
      CACHE INTERNAL "TCL include directories"
      FORCE)
ELSE(UNIX OR APPLE)
  IF(WIN32)
    SET(SCI_TCL_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/generic;${CMAKE_CURRENT_SOURCE_DIR}/win"
        CACHE INTERNAL "TCL include directories"
        FORCE)
  ENDIF(WIN32)
ENDIF(UNIX)

SET(SCI_TCLSH tclsh CACHE INTERNAL "Tclsh executable" FORCE)
SET(TCLPIPEDLLNAME "tclpip${TCL_MAJOR_VERSION}${TCL_MINOR_VERSION}"
    CACHE INTERNAL "Tcl pipe library name"
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
ELSE(WIN32 AND MSVC)
  IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    STRING(REGEX REPLACE "-O3" "-O0" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
  ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
ENDIF(WIN32 AND MSVC)

#################################################################
# Building library

ADD_LIBRARY(${SCI_TCLSTUB_LIBRARY} STATIC ${TCL_STUB_SRCS})

SET_TARGET_PROPERTIES(${SCI_TCLSTUB_LIBRARY}
                      PROPERTIES
                        OUTPUT_NAME tclstub.${TCL_VERSION}
                        COMPILE_DEFINITIONS "STATIC_BUILD"
                        COMPILE_DEFINITIONS_Debug "TCL_CFG_DEBUG=1"
                     )
IF(UNIX AND NOT APPLE AND NOT CYGWIN)
  SET_TARGET_PROPERTIES(${SCI_TCLSTUB_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fPIC"
                       )
ENDIF(UNIX AND NOT APPLE AND NOT CYGWIN)

ADD_LIBRARY(${SCI_TCL_LIBRARY} ${TCL_SRCS})
SET_TARGET_PROPERTIES(${SCI_TCL_LIBRARY}
                      PROPERTIES
                        OUTPUT_NAME tcl.${TCL_VERSION}
                        COMPILE_DEFINITIONS_DEBUG "TCL_CFG_DEBUG=1"
                     )
IF(WIN32 AND MSVC)
  # Turn off Visual Studio warnings about strcpy et.al. being unsafe
  # for release builds (only supporting CMake defaults).
  SET_TARGET_PROPERTIES(${SCI_TCL_LIBRARY}
                        PROPERTIES
                          COMPILE_DEFINITIONS_RELEASE "_CRT_SECURE_NO_WARNINGS=1"
                          COMPILE_DEFINITIONS_MINSIZEREL "_CRT_SECURE_NO_WARNINGS=1"
                       )

ENDIF(WIN32 AND MSVC)

IF(ENABLE_PIC)
  SET_TARGET_PROPERTIES(${SCI_TCL_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fPIC")
ENDIF(ENABLE_PIC)

IF(WIN32)
  ADD_DEFINITIONS(-DTCL_PIPE_DLL=\"${TCLPIPEDLLNAME}\")

  # Windows socket 2 library
  FIND_LIBRARY(WS2_LIBRARY ws2_32 PATHS ${MS_SDK_REG_PATH})

  ADD_LIBRARY(${TCLPIPEDLLNAME} ${TCL_WIN32_PIPE_SRCS})
  TARGET_LINK_LIBRARIES(${TCLPIPEDLLNAME}
                        ${WS2_LIBRARY}
                       )

  TARGET_LINK_LIBRARIES(${SCI_TCL_LIBRARY}
                        ${SCI_ZLIB_LIBRARY}
                        ${WS2_LIBRARY}
                       )
  SET_TARGET_PROPERTIES(${SCI_TCL_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fp:strict -Dinline=__inline"
                          COMPILE_DEFINITIONS "TCL_THREADS;BUILD_tcl;HAVE_ZLIB=1;_MBCS"
                          DEFINE_SYMBOL _USRDLL
                          DEFINE_SYMBOL _WINDLL
                       )
ELSE(WIN32)
  
  FIND_PACKAGE( Threads )
  SET(TCL_THREAD_LIBRARY ${CMAKE_THREAD_LIBS})          
  ADD_DEFINITIONS(-DTCL_THREADS)
  ADD_DEFINITIONS(-DHAVE_ZLIB=1)
  TARGET_LINK_LIBRARIES(${SCI_TCL_LIBRARY}
                        ${SCI_ZLIB_LIBRARY}
                        ${TCL_THREAD_LIBRARY}
                       )  
ENDIF(WIN32)


IF(UNIX)
  # These libraries are not needed under Windows.
  FIND_LIBRARY(M_LIBRARY m)
  FIND_LIBRARY(DL_LIBRARY dl)
ENDIF(UNIX)


ADD_EXECUTABLE(${SCI_TCLSH} ${TCLSH_SRCS})
IF(UNIX)
  TARGET_LINK_LIBRARIES(${SCI_TCLSH} ${SCI_TCL_LIBRARY} ${M_LIBRARY} ${DL_LIBRARY})
ELSE(UNIX)
  IF(WIN32)
    TARGET_LINK_LIBRARIES(${SCI_TCLSH} ${SCI_TCL_LIBRARY})
  ENDIF(WIN32)
ENDIF(UNIX)

IF(APPLE)
  SET_TARGET_PROPERTIES(${SCI_TCL_LIBRARY} PROPERTIES LINK_FLAGS  "-framework CoreFoundation")
  SET_TARGET_PROPERTIES(${SCI_TCLSH} PROPERTIES LINK_FLAGS  "-framework CoreFoundation")
ENDIF(APPLE)
