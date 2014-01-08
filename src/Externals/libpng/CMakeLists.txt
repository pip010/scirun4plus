# libpng source code license
#
# This copy of the libpng notices is provided for your convenience.  In case of
# any discrepancy between this copy and the notices in the file png.h that is
# included in the libpng distribution, the latter shall prevail.
# 
# COPYRIGHT NOTICE, DISCLAIMER, and LICENSE:
# 
# If you modify libpng you may insert additional notices immediately following
# this sentence.
# 
# libpng versions 1.2.6, August 15, 2004, through 1.2.35, February 14, 2009, are
# Copyright (c) 2004, 2006-2008 Glenn Randers-Pehrson, and are
# distributed according to the same disclaimer and license as libpng-1.2.5
# with the following individual added to the list of Contributing Authors
# 
#    Cosmin Truta
# 
# libpng versions 1.0.7, July 1, 2000, through 1.2.5 - October 3, 2002, are
# Copyright (c) 2000-2002 Glenn Randers-Pehrson, and are
# distributed according to the same disclaimer and license as libpng-1.0.6
# with the following individuals added to the list of Contributing Authors
# 
#    Simon-Pierre Cadieux
#    Eric S. Raymond
#    Gilles Vollant
# 
# and with the following additions to the disclaimer:
# 
#    There is no warranty against interference with your enjoyment of the
#    library or against infringement.  There is no warranty that our
#    efforts or the library will fulfill any of your particular purposes
#    or needs.  This library is provided with all faults, and the entire
#    risk of satisfactory quality, performance, accuracy, and effort is with
#    the user.
# 
# libpng versions 0.97, January 1998, through 1.0.6, March 20, 2000, are
# Copyright (c) 1998, 1999 Glenn Randers-Pehrson, and are
# distributed according to the same disclaimer and license as libpng-0.96,
# with the following individuals added to the list of Contributing Authors:
# 
#    Tom Lane
#    Glenn Randers-Pehrson
#    Willem van Schaik
# 
# libpng versions 0.89, June 1996, through 0.96, May 1997, are
# Copyright (c) 1996, 1997 Andreas Dilger
# Distributed according to the same disclaimer and license as libpng-0.88,
# with the following individuals added to the list of Contributing Authors:
# 
#    John Bowler
#    Kevin Bracey
#    Sam Bushell
#    Magnus Holmgren
#    Greg Roelofs
#    Tom Tanner
# 
# libpng versions 0.5, May 1995, through 0.88, January 1996, are
# Copyright (c) 1995, 1996 Guy Eric Schalnat, Group 42, Inc.
# 
# For the purposes of this copyright and license, "Contributing Authors"
# is defined as the following set of individuals:
# 
#    Andreas Dilger
#    Dave Martindale
#    Guy Eric Schalnat
#    Paul Schmidt
#    Tim Wegner
# 
# The PNG Reference Library is supplied "AS IS".  The Contributing Authors
# and Group 42, Inc. disclaim all warranties, expressed or implied,
# including, without limitation, the warranties of merchantability and of
# fitness for any purpose.  The Contributing Authors and Group 42, Inc.
# assume no liability for direct, indirect, incidental, special, exemplary,
# or consequential damages, which may result from the use of the PNG
# Reference Library, even if advised of the possibility of such damage.
# 
# Permission is hereby granted to use, copy, modify, and distribute this
# source code, or portions hereof, for any purpose, without fee, subject
# to the following restrictions:
# 
# 1. The origin of this source code must not be misrepresented.
# 
# 2. Altered versions must be plainly marked as such and must not
#    be misrepresented as being the original source.
# 
# 3. This Copyright notice may not be removed or altered from any
#    source or altered source distribution.
# 
# The Contributing Authors and Group 42, Inc. specifically permit, without
# fee, and encourage the use of this source code as a component to
# supporting the PNG file format in commercial products.  If you use this
# source code in a product, acknowledgment is not required but would be
# appreciated.
# 
# 
# A "png_get_copyright" function is available, for convenient use in "about"
# boxes and the like:
# 
#    printf("%s",png_get_copyright(NULL));
# 
# Also, the PNG logo (in PNG format, of course) is supplied in the
# files "pngbar.png" and "pngbar.jpg (88x31) and "pngnow.png" (98x31).
# 
# Libpng is OSI Certified Open Source Software.  OSI Certified Open Source is a
# certification mark of the Open Source Initiative.
# 
# Glenn Randers-Pehrson
# glennrp at users.sourceforge.net
# February 14, 2009

# CMake build scripts modifications by Scientific Computing and Imaging Institute
# Copyright (c) 2009 Scientific Computing and Imaging Institute, University of Utah.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(VERSION 2.6)
ENDIF(COMMAND CMAKE_POLICY)

project(PNG)

IF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)
  SET(ENABLE_PIC ON)
ENDIF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)

# Copyright (C) 2007 Glenn Randers-Pehrson
# For conditions of distribution and use, see copyright notice in png.h

set(PNGLIB_MAJOR 1)
set(PNGLIB_MINOR 2)
set(PNGLIB_RELEASE 35)
set(PNGLIB_NAME libpng${PNGLIB_MAJOR}${PNGLIB_MINOR})
set(PNGLIB_VERSION ${PNGLIB_MAJOR}.${PNGLIB_MINOR}.${PNGLIB_RELEASE})

# needed packages

## Try out using a path set in the SCIRun CMakeList.txt file
if(DEFINED SCI_ZLIB_LIBRARY)

  SET(ZLIB_LIBRARY ${SCI_ZLIB_LIBRARY})
  SET(ZLIB_INCLUDE ${SCI_ZLIB_INCLUDE})

else(DEFINED SCI_ZLIB_LIBRARY)

  # try to find a system zlib
  find_package(ZLIB REQUIRED)

endif(DEFINED SCI_ZLIB_LIBRARY)

if(NOT WIN32)
 find_library(M_LIBRARY
     NAMES m
     PATHS /usr/lib /usr/local/lib
 )
 if(NOT M_LIBRARY)
   message(STATUS
     "math library 'libm' not found - floating point support disabled")
 endif(NOT M_LIBRARY)
else(NOT WIN32)
 # not needed on windows
 set(M_LIBRARY "")
endif(NOT WIN32)


# TODO: check this
# maybe needs improving, but currently I don't know when we can enable what :)
set(png_asm_tmp "OFF")
if(NOT WIN32)
 find_program(uname_executable NAMES uname PATHS /bin /usr/bin /usr/local/bin)
 if(uname_executable)
   EXEC_PROGRAM(${uname_executable} ARGS --machine OUTPUT_VARIABLE uname_output)
   if("uname_output" MATCHES "^.*i[1-9]86.*$")
      set(png_asm_tmp "ON")
   else("uname_output" MATCHES "^.*i[1-9]86.*$")
      set(png_asm_tmp "OFF")
   endif("uname_output" MATCHES "^.*i[1-9]86.*$")
 endif(uname_executable)
else(NOT WIN32)
 # this env var is normally only set on win64
 SET(TEXT "ProgramFiles(x86)")
 if("$ENV{${TEXT}}" STREQUAL "")
  set(png_asm_tmp "ON")
 endif("$ENV{${TEXT}}" STREQUAL "")
endif(NOT WIN32)

SET(SCI_PNG_LIBRARY png CACHE INTERNAL "PNG library name." FORCE)
SET(SCI_PNG_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "Public libpng include directories." FORCE)

# OUR SOURCES
set(libpng_sources
 png.h
 pngconf.h
 png.c
 pngerror.c
 pngget.c
 pngmem.c
 pngpread.c
 pngread.c
 pngrio.c
 pngrtran.c
 pngrutil.c
 pngset.c
 pngtrans.c
 pngwio.c
 pngwrite.c
 pngwtran.c
 pngwutil.c
)

# SOME NEEDED DEFINITIONS
if(WIN32 AND MSVC)
  add_definitions(-DPNG_NO_MODULEDEF -D_BIND_TO_CURRENT_VCLIBS_VERSION=1)

  # MSVC level 3 warnings are a little ridiculous
  SET(CMAKE_C_WARNING_LEVEL 1)
  IF(CMAKE_C_FLAGS MATCHES "/W[0-4]")
    STRING(REGEX REPLACE "/W[0-4]" "/W1"
            CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  ELSE(CMAKE_C_FLAGS MATCHES "/W[0-4]")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W1")
  ENDIF(CMAKE_C_FLAGS MATCHES "/W[0-4]")
endif(WIN32 AND MSVC)

add_definitions(-DZLIB_DLL)

add_definitions(-DLIBPNG_NO_MMX)
add_definitions(-DPNG_NO_MMX_CODE)
add_definitions(-DPNG_CONSOLE_IO_SUPPORTED)

if(NOT M_LIBRARY AND NOT WIN32)
 add_definitions(-DPNG_NO_FLOATING_POINT_SUPPORTED)
endif(NOT M_LIBRARY AND NOT WIN32)

# NOW BUILD OUR TARGET
include_directories(${PNG_SOURCE_DIR} ${ZLIB_INCLUDE})

add_library(${SCI_PNG_LIBRARY} STATIC ${libpng_sources})
set_target_properties(${SCI_PNG_LIBRARY}
                      PROPERTIES
                        OUTPUT_NAME png${PNGLIB_MAJOR}${PNGLIB_MINOR}
                     )
target_link_libraries(${SCI_PNG_LIBRARY} ${ZLIB_LIBRARY} ${M_LIBRARY})

if(WIN32) 
  set_target_properties(${SCI_PNG_LIBRARY}
                        PROPERTIES
                          DEFINE_SYMBOL PNG_BUILD_DLL
                       )
endif(WIN32) 

if(ENABLE_PIC)
  set_target_properties(${SCI_PNG_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fPIC"
                       )
endif(ENABLE_PIC)
