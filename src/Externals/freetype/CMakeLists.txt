# FreeType 2 configuration rules for UNIX platforms
#
# Copyright 1996-2000, 2002, 2004, 2006 by
# David Turner, Robert Wilhelm, and Werner Lemberg.
#
# This file is part of the FreeType project, and may only be used, modified,
# and distributed under the terms of the FreeType project license,
# LICENSE.TXT.  By continuing to use, modify, or distribute this file you
# indicate that you have read the license and understand and accept it
# fully.

# CMake build scripts provided by Scientific Computing and Imaging Institute
# Copyright (c) 2009 Scientific Computing and Imaging Institute, University of Utah.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(VERSION 2.6)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(freetype)
#-----------------------------------------------------------------------------

IF(BUILD_STATIC_THIRDPARTY)
  IF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)
    SET(ENABLE_PIC ON)
  ENDIF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)
  SET(BUILD_SHARED_LIBS OFF)
ENDIF(BUILD_STATIC_THIRDPARTY)

# Test for some required system information.

INCLUDE (CMakeBackwardCompatibilityC)

# will have to change freetype config headers to export DLL entry point symbol

SET(FREETYPE_autofit_SRCS
     src/autofit/autofit.c
   )

SET(FREETYPE_base_SRCS
     src/base/ftbase.c
     src/base/ftbbox.c
     src/base/ftbdf.c
     src/base/ftbitmap.c
     src/base/ftdebug.c
     src/base/ftglyph.c
     src/base/ftinit.c
     src/base/ftmm.c
     src/base/ftotval.c
     src/base/ftpfr.c
     src/base/ftstroke.c
     src/base/ftsynth.c
     src/base/ftsystem.c
     src/base/fttype1.c
     src/base/ftwinfnt.c
     src/base/ftxf86.c
   )

SET(FREETYPE_bdf_SRCS
     src/bdf/bdf.c
   )

SET(FREETYPE_cache_SRCS
     src/cache/ftcache.c
   )

SET(FREETYPE_cff_SRCS
     src/cff/cff.c
   )

SET(FREETYPE_cid_SRCS
     src/cid/type1cid.c
   )

SET(FREETYPE_ftgzip_SRCS
     src/gzip/ftgzip.c
   )

SET(FREETYPE_lzw_SRCS
     src/lzw/ftlzw.c
   )

SET(FREETYPE_otvalid_SRCS
     src/otvalid/otvalid.c
   )

SET(FREETYPE_pcf_SRCS
     src/pcf/pcf.c
   )

SET(FREETYPE_pfr_SRCS
     src/pfr/pfr.c
   )

SET(FREETYPE_psaux_SRCS
     src/psaux/psaux.c
   )

SET(FREETYPE_pshinter_SRCS
     src/pshinter/pshinter.c
   )

SET(FREETYPE_psnames_SRCS
     src/psnames/psnames.c
   )

SET(FREETYPE_raster_SRCS
     src/raster/raster.c
   )

SET(FREETYPE_sfnt_SRCS
     src/sfnt/sfnt.c
   )

SET(FREETYPE_smooth_SRCS
     src/smooth/smooth.c
   )

SET(FREETYPE_truetype_SRCS
     src/truetype/truetype.c
   )

SET(FREETYPE_type1_SRCS
     src/type1/type1.c
   )

SET(FREETYPE_type42_SRCS
     src/type42/type42.c
   )

SET(FREETYPE_winfonts_SRCS
     src/winfonts/winfnt.c
   )

SET(FREETYPE_SRCS
     ${FREETYPE_autofit_SRCS}
     ${FREETYPE_base_SRCS}
     ${FREETYPE_bdf_SRCS}
     ${FREETYPE_cache_SRCS}
     ${FREETYPE_cff_SRCS}
     ${FREETYPE_cid_SRCS}
     ${FREETYPE_ftgzip_SRCS}
     ${FREETYPE_lzw_SRCS}
     ${FREETYPE_otvalid_SRCS}
     ${FREETYPE_pcf_SRCS}
     ${FREETYPE_pfr_SRCS}
     ${FREETYPE_psaux_SRCS}
     ${FREETYPE_pshinter_SRCS}
     ${FREETYPE_psnames_SRCS}
     ${FREETYPE_raster_SRCS}
     ${FREETYPE_sfnt_SRCS}
     ${FREETYPE_smooth_SRCS}
     ${FREETYPE_truetype_SRCS}
     ${FREETYPE_type1_SRCS}
     ${FREETYPE_type42_SRCS}
     ${FREETYPE_winfonts_SRCS}
   )

SET(SCI_FREETYPE_INCLUDE
      ${CMAKE_CURRENT_SOURCE_DIR}/include
    CACHE INTERNAL
    "Public freetype include directories."
    FORCE
   )

SET(SCI_FREETYPE_LIBRARY freetype CACHE INTERNAL "Freetype library name." FORCE)

ADD_DEFINITIONS(-DFT2_BUILD_LIBRARY)
IF(APPLE)
  ADD_DEFINITIONS(-DDARWIN_NO_CARBON)
ENDIF(APPLE)

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

INCLUDE_DIRECTORIES(${SCI_FREETYPE_INCLUDE} ${SCI_ZLIB_INCLUDE})
ADD_LIBRARY(${SCI_FREETYPE_LIBRARY} ${FREETYPE_SRCS}) 
IF(WIN32)
  SET_TARGET_PROPERTIES(${SCI_FREETYPE_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-DWIN32 -D_LIB -D_MBCS"
                          DEFINE_SYMBOL BUILD_DLL)
ENDIF(WIN32)

IF(ENABLE_PIC)
  SET_TARGET_PROPERTIES(${SCI_FREETYPE_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fPIC")
ENDIF(ENABLE_PIC)


TARGET_LINK_LIBRARIES(${SCI_FREETYPE_LIBRARY} ${SCI_ZLIB_LIBRARY})
