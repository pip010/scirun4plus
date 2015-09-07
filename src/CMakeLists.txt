#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2010 Scientific Computing and Imaging Institute,
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

# -- Determine the kernel version of OSX
# -- 8 and less are OSX 10.0 - 10.4 (no longer officially supported)
# -- 9.x is 10.5 (Leopard)
# -- 10.x is 10.6 (Snow Leopard)
# -- 11.x is 10.7 (Lion)
#
# TODO: this matches ALL OS X versions > kernel 8.x
# Is it still needed? If so, should it be used with OS X 10.6 and higher? 
# For now, assume APPLE_LEOPARD refers any version of OS X from Leopard onwards
IF (APPLE)
  EXEC_PROGRAM(uname ARGS -v  OUTPUT_VARIABLE DARWIN_VERSION)
  STRING(REGEX MATCH "[0-9]+" DARWIN_MAJOR_VERSION ${DARWIN_VERSION})
  MESSAGE(STATUS "Detecting Darwin kernel: DARWIN_MAJOR_VERSION=${DARWIN_MAJOR_VERSION}")
  IF (DARWIN_MAJOR_VERSION GREATER 8)
    SET(APPLE_LEOPARD 1 INTERNAL)
    ADD_DEFINITIONS(-DAPPLE_LEOPARD)
    IF (DARWIN_MAJOR_VERSION LESS 11)
      SET(APPLE_LEOPARD_OR_SNOW_LEOPARD 1 INTERNAL)
      ADD_DEFINITIONS(-DAPPLE_LEOPARD_OR_SNOW_LEOPARD)
    ENDIF (DARWIN_MAJOR_VERSION LESS 11)
  ELSEIF (DARWIN_MAJOR_VERSION GREATER 8)
    MESSAGE(WARNING "SCIRun supports Mac OS X 10.5 and higher, and is not guaranteed to work on older versions of the operating system.")
  ENDIF (DARWIN_MAJOR_VERSION GREATER 8)
ENDIF(APPLE)


IF(WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
  ADD_DEFINITIONS(-D_BIND_TO_CURRENT_VCLIBS_VERSION=1)
ENDIF(WIN32)

###############################################################################
# Don't allow the user to build in the src directory...
#

IF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  MESSAGE( "" )
  MESSAGE( "Building SCIRun in the source directory is not allowed." )
  MESSAGE( "(In order to clean up, you now need to run: rm -rf CMakeFiles CMakeCache.txt" )
  MESSAGE( "" )
  MESSAGE( FATAL_ERROR "Goodbye." )
ENDIF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})

CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(VERSION 2.8)
ENDIF(COMMAND CMAKE_POLICY)


PROJECT(SCIRUN_CORE)

# SOME HOW CMAKE DOES NOT GET IT RIGHT
IF(APPLE)
IF(CMAKE_OSX_ARCHITECTURES MATCHES .*64)
  SET(CMAKE_SIZEOF_VOID_P 8)
ENDIF(CMAKE_OSX_ARCHITECTURES MATCHES .*64)
ENDIF(APPLE)

########################################################################
# Test for some required system information.
INCLUDE (CMakeBackwardCompatibilityC)

########################################################################
# SCIRUN version number.

SET(SCIRUN_VERSION_MAJOR "4")
SET(SCIRUN_VERSION_MINOR "7")
SET(SCIRUN_VERSION_PATCH "0")

SET(SCIRUN_VERSION_STRING "${SCIRUN_VERSION_MAJOR}.${SCIRUN_VERSION_MINOR}")

########################################################################
# Configure testing support.  This should be done before any
# MESSAGE(FATAL_ERROR ...) commands are invoked.

INCLUDE( CTest )

IF(BUILD_TESTING)
  ENABLE_TESTING()
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)

SET(SCIRUN_CORE_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

########################################################################
# Output directories.

SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${SCIRUN_CORE_BINARY_DIR}/lib CACHE INTERNAL "Single output directory for building all libraries.")
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SCIRUN_CORE_BINARY_DIR} CACHE INTERNAL "Single output directory for building all executables.")
MARK_AS_ADVANCED(CMAKE_LIBRARY_OUTPUT_DIRECTORY CMAKE_RUNTIME_OUTPUT_DIRECTORY)

SET(BUILD_FROM_TOP_LEVEL 1)

IF(NOT WIN32 AND NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT WIN32 AND NOT CMAKE_BUILD_TYPE)

SET(CMAKE_VERBOSE_MAKEFILE off CACHE STRING "Verbose makefile output")
MARK_AS_ADVANCED(CMAKE_VERBOSE_MAKEFILE)

########################################################################
# shared/static

OPTION(BUILD_SHARED_LIBS "Build shared libraries" ON)
MARK_AS_ADVANCED(BUILD_SHARED_LIBS)

OPTION(BUILD_STATIC_THIRDPARTY "Build third party libraries static" OFF)

########################################################################
# Build SCIRun dataflow, executable and packages by default.
# If this is not set we just build the core libraries.
# Primarily used by Seg3D build script (deprecated).

OPTION(BUILD_DATAFLOW "Build SCIRun executable, dataflow and packages" ON)
MARK_AS_ADVANCED(BUILD_DATAFLOW)

########################################################################
# Build BioMesh3D pipeline.
OPTION(BUILD_BIOMESH3D "Build BioMesh3D pipeline" ON)

 IF(BUILD_BIOMESH3D)
   ADD_SUBDIRECTORY(Applications/FEMesher)
 ENDIF(BUILD_BIOMESH3D)


########################################################################
# Build BioMesh3D Remote Client.
OPTION(BUILD_BIOMESH3D_REMOTE_SUPPORT "Build BioMesh3D Remote Support." OFF)

 IF(BUILD_BIOMESH3D_REMOTE_SUPPORT)
   ADD_DEFINITIONS(-DBUILD_BIOMESH3D_REMOTE_SUPPORT)
 ENDIF(BUILD_BIOMESH3D_REMOTE_SUPPORT)

########################################################################
#  Tetgen may also be used outside of the meshing pipeline.

OPTION(WITH_TETGEN "Build Tetgen." ON)
MARK_AS_ADVANCED(WITH_TETGEN)
IF(BUILD_BIOMESH3D)
  SET(WITH_TETGEN ON)
ENDIF(BUILD_BIOMESH3D)

########################################################################
#  Build Insight Toolkit.  Turned off by default, must be on for Seg3D.

OPTION(WITH_ITK "Build with Insight Toolkit" OFF)
MARK_AS_ADVANCED(WITH_ITK)		   

########################################################################
# GDCM Library (which could be pulled from ITK or specified individually)
OPTION(WITH_GDCM "Compile with GDCM support" OFF)

IF(WITH_ITK OR GDCM_DIR)
  SET(WITH_GDCM ON)
ENDIF(WITH_ITK OR GDCM_DIR)
MARK_AS_ADVANCED(FORCE WITH_GDCM)

# use provided GDCM Library instead of ITK's version
IF(GDCM_DIR)
  SET(ITK_USE_SYSTEM_GDCM ON)
ELSE(GDCM_DIR)
  SET(ITK_USE_SYSTEM_GDCM OFF)
ENDIF(GDCM_DIR)

########################################################################
# LAPACK and BLAS libraries
# TODO: replace with customized CMake module
#
# search for LAPACK

OPTION(WITH_LAPACK "Build with LAPACK library (Note: we use the vecLib framework on OS X by default)." ON)
MARK_AS_ADVANCED(CLEAR WITH_LAPACK)

OPTION(BUILD_NETLIB "Build SCIRun's default LAPACK library. Do NOT enable if you want to use external LAPACK libraries." OFF)
MARK_AS_ADVANCED(BUILD_NETLIB)

IF(WITH_LAPACK)
  IF(APPLE)
    # OS X has one library (vecLib) for both BLAS and LAPACK
    SET(WITH_BLAS ON)
    SET(HAVE_LAPACK "YES")
    SET(LAPACK_LIBRARY "-framework vecLib" CACHE STRING "Lapack library name" FORCE)
  ELSE(APPLE)
    # everybody else
    IF(BUILD_NETLIB)
      # Build LAPACK distribution in src/Externals,
      MESSAGE(STATUS "Building default SCIRun LAPACK library.")
      SET(HAVE_LAPACK "YES")

    ELSE(BUILD_NETLIB)
      MESSAGE(STATUS "Looking for LAPACK")
      SET(LAPACK_DIR "" CACHE PATH "Path to LAPACK" FORCE)
      MARK_AS_ADVANCED(CLEAR LAPACK_DIR)

      FIND_LIBRARY(LAPACK_LIBRARY NAMES lapack clapack flapack
                   PATHS ${LAPACK_DIR}
                     ${LAPACK_DIR}/lib NO_DEFAULT_PATH)
      FIND_LIBRARY(LAPACK_LIBRARY NAMES lapack clapack flapack
                   PATHS
                     /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64 /opt/lib /opt/lib64 /lib  /lib64 ENV LD_LIBRARY_PATH)

      IF(LAPACK_LIBRARY)
        MESSAGE(STATUS "Looking for LAPACK - found (${LAPACK_LIBRARY})")
        SET(HAVE_LAPACK "YES")

      ELSE(LAPACK_LIBRARY)
        SET(BUILD_NETLIB ON CACHE BOOL "Build SCIRun's default LAPACK and BLAS libraries. Do NOT enable if you want to use external LAPACK libraries." FORCE)
        MESSAGE(STATUS "Looking for LAPACK - not found. Building default SCIRun LAPACK library.")
        SET(HAVE_LAPACK "YES")
      ENDIF(LAPACK_LIBRARY)
    ENDIF(BUILD_NETLIB)
  ENDIF(APPLE)
ELSE(WITH_LAPACK)
  SET(HAVE_LAPACK "NO")
  UNSET(LAPACK_LIBRARY CACHE)
ENDIF(WITH_LAPACK)

########################################################################
# search for BLAS

OPTION(WITH_BLAS "Build with BLAS library" OFF)
MARK_AS_ADVANCED(WITH_BLAS)

IF(WITH_BLAS)
  MARK_AS_ADVANCED(CLEAR WITH_BLAS)
ENDIF(WITH_BLAS)

IF(WITH_BLAS)
  IF(APPLE)
    # TODO: setting atlas library may be redundant    
    # SET(ATLAS_LIBRARY "-framework vecLib" CACHE STRING "Atlas library name" FORCE)

    # OS X has one library (vecLib) for both BLAS and LAPACK:
    SET(BLAS_LIBRARIES "-framework vecLib" CACHE STRING "Blas library name" FORCE)
  ELSE(APPLE)
    MESSAGE(STATUS "Looking for BLAS")
    SET(BLAS_DIR "" CACHE PATH "Path to BLAS" FORCE)
    MARK_AS_ADVANCED(CLEAR BLAS_DIR)

    IF(BLAS_DIR)
      FIND_LIBRARY(CBLAS_LIBRARY NAMES cblas gslcblas PATHS ${BLAS_DIR} ${BLAS_DIR}/lib)
      FIND_LIBRARY(FBLAS_LIBRARY NAMES fblas f77blas blas PATHS ${BLAS_DIR} ${BLAS_DIR}/lib)
      IF(NOT CBLAS_LIBRARY)
        # ATLAS provides a complete BLAS API (C, Fortran) and a very small subset of the LAPACK API
        FIND_LIBRARY(ATLAS_LIBRARY NAMES atlas PATHS ${BLAS_DIR} ${BLAS_DIR}/lib)
      ENDIF(NOT CBLAS_LIBRARY)

      FIND_PATH(CBLAS_INCLUDE cblas.h ${BLAS_DIR} ${BLAS_DIR}/include ${BLAS_DIR}/../include)
    ELSE(BLAS_DIR)
      IF(LAPACK_DIR)
        # BLAS and LAPACK are traditionally in the same place
        FIND_LIBRARY(CBLAS_LIBRARY NAMES cblas gslcblas PATHS ${LAPACK_DIR} ${LAPACK_DIR}/lib)
        FIND_LIBRARY(FBLAS_LIBRARY NAMES fblas f77blas blas PATHS ${LAPACK_DIR} ${LAPACK_DIR}/lib)

        IF(NOT CBLAS_LIBRARY)
          FIND_LIBRARY(ATLAS_LIBRARY NAMES atlas PATHS ${LAPACK_DIR} ${LAPACK_DIR}/lib)
        ENDIF(NOT CBLAS_LIBRARY)

        FIND_PATH(CBLAS_INCLUDE cblas.h ${LAPACK_DIR} ${LAPACK_DIR}/include)

      ELSE(LAPACK_DIR)

        SET(SEARCH_PATHS /usr/local/lib /usr/lib /opt/lib /lib /usr/local/lib64 /usr/lib64 /opt/lib64 /lib64)
        FIND_LIBRARY(CBLAS_LIBRARY NAMES cblas gslcblas PATHS ${SEARCH_PATHS} ENV LD_LIBRARY_PATH)
        FIND_LIBRARY(FBLAS_LIBRARY NAMES fblas f77blas blas PATHS ${SEARCH_PATHS} ENV LD_LIBRARY_PATH)
        IF(NOT CBLAS_LIBRARY)
          FIND_LIBRARY(ATLAS_LIBRARY NAMES atlas PATHS ${SEARCH_PATHS} ENV LD_LIBRARY_PATH)
        ENDIF(NOT CBLAS_LIBRARY)

        FIND_PATH(CBLAS_INCLUDE cblas.h /usr/local/include/ /usr/include /opt/include /include)
      ENDIF(LAPACK_DIR)
    ENDIF(BLAS_DIR)

    IF(CBLAS_LIBRARY AND CBLAS_INCLUDE)
      MESSAGE(STATUS "Looking for BLAS - found CBLAS")
      SET(HAVE_CBLAS "YES" CACHE INTERNAL "Do not edit" FORCE)

      SET(BLAS_LIBRARIES ${CBLAS_LIBRARY})
      SET(SCIRUN_INCLUDE_THIRDPARTY ${SCIRUN_INCLUDE_THIRDPARTY} ${BLAS_INCLUDE})
    ELSE(CBLAS_LIBRARY AND CBLAS_INCLUDE)
      IF(ATLAS_LIBRARY AND CBLAS_INCLUDE)
        MESSAGE(STATUS "Looking for BLAS - found ATLAS")
        SET(HAVE_BLAS "YES" CACHE INTERNAL "Do not edit" FORCE)
        SET(HAVE_CBLAS "YES" CACHE INTERNAL "Do not edit" FORCE)

        SET(BLAS_LIBRARIES ${ATLAS_LIBRARY})
        SET(SCIRUN_INCLUDE_THIRDPARTY ${SCIRUN_INCLUDE_THIRDPARTY} ${BLAS_INCLUDE})
      ELSE(ATLAS_LIBRARY AND CBLAS_INCLUDE)
        MESSAGE(STATUS "BLAS and ATLAS not found - trying FBLAS")
      ENDIF(ATLAS_LIBRARY AND CBLAS_INCLUDE)
    ENDIF(CBLAS_LIBRARY AND CBLAS_INCLUDE)

    IF(FBLAS_LIBRARY)
      MESSAGE(STATUS "Looking for BLAS - found FBLAS")
      SET(HAVE_BLAS "YES" CACHE INTERNAL "Do not edit" FORCE)
      LIST(APPEND BLAS_LIBRARIES ${FBLAS_LIBRARY})
    ELSE(FBLAS_LIBRARY)
      MESSAGE(STATUS "Looking for BLAS - did not find FBLAS")
    ENDIF(FBLAS_LIBRARY)

    MARK_AS_ADVANCED(BLAS_LIBRARIES)
    MARK_AS_ADVANCED(CBLAS_LIBRARY)
    MARK_AS_ADVANCED(FBLAS_LIBRARY)
    MARK_AS_ADVANCED(ATLAS_LIBRARY)
    MARK_AS_ADVANCED(CBLAS_INCLUDE)
  ENDIF(APPLE)
ELSE(WITH_BLAS)
  SET(HAVE_CBLAS "NO")
  SET(WITH_BLAS OFF)
  UNSET(BLAS_LIBRARIES CACHE)
ENDIF(WITH_BLAS)

########################################################################
# Compiler flags

IF(CMAKE_COMPILER_IS_GNUCC)
  # Always turn on more verbose warnings when compiling with gcc and g++.
  SET(CMAKE_CXX_FLAGS "-Wall -fno-strict-aliasing ${CMAKE_CXX_FLAGS}")
  SET(CMAKE_C_FLAGS "-Wall -fno-strict-aliasing ${CMAKE_C_FLAGS}")

  # Additional warnings for debug builds
  SET(CMAKE_CXX_FLAGS_DEBUG "-Wold-style-cast ${CMAKE_CXX_FLAGS_DEBUG}")
  IF(APPLE)
    SET(CMAKE_CXX_FLAGS_DEBUG "-Wshorten-64-to-32 ${CMAKE_CXX_FLAGS_DEBUG}")
  ENDIF(APPLE)
ENDIF(CMAKE_COMPILER_IS_GNUCC)

IF(WIN32 AND MSVC)
  SET(CMAKE_CXX_FLAGS "/MP ${CMAKE_CXX_FLAGS}") 
ENDIF(WIN32 AND MSVC)

################################################################################################
# Swichting off OpenGL support
#

SET(DISABLE_OPENGL OFF CACHE BOOL "Switch off OpenGL Support in SCIRun")
MARK_AS_ADVANCED (DISABLE_OPENGL)

IF (NOT DISABLE_OPENGL)
  SET(DEF_OPENGL "#define HAVE_OPENGL 1" CACHE INTERNAL "internal variable" FORCE)
ELSE (NOT DISABLE_OPENGL)
  SET(DEF_OPENGL "" CACHE INTERNAL "internal variable" FORCE)
ENDIF (NOT DISABLE_OPENGL)

################################################################################################
# Number of Bits (32 | 64) Check
#

IF(CMAKE_SIZEOF_VOID_P MATCHES  8)
  SET(NBITS 64)
  ADD_DEFINITIONS(-DSCI_64BITS)
ELSE(CMAKE_SIZEOF_VOID_P MATCHES 8)
  SET(NBITS 32)
  ADD_DEFINITIONS(-DSCI_32BITS)
ENDIF(CMAKE_SIZEOF_VOID_P MATCHES 8)

MESSAGE(STATUS "FYI: This is a ${NBITS} bit build.")

#MARK_AS_ADVANCED(CLEAR SCIRUN_ENABLE_64BIT)

################################################################################################
# Under OSX, use single_module (as opposed to multiple_module) linking.

if(APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "-single_module" CACHE STRING "Flags used to link a shared library." FORCE)
endif(APPLE)

#############################################################################################
# Installation helpers

# Run this on executables you wish to be installed
MACRO(SCIRUN_INSTALL_EXECUTABLE executable)
  #   SET_TARGET_PROPERTIES(${executable} PROPERTIES
  #     INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/bin
  #     )
  INSTALL(TARGETS ${executable}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
ENDMACRO(SCIRUN_INSTALL_EXECUTABLE)

MACRO(SCIRUN_INSTALL_LIBRARY library)
  #   IF(APPLE)
  #     SET_TARGET_PROPERTIES(${library} PROPERTIES
  #       BUILD_WITH_INSTALL_RPATH ON
  #       INSTALL_NAME_DIR "@executable_path/../lib"
  #       )
  #   ENDIF(APPLE)

  # TODO(bigler): don't install static libraries?
  INSTALL(TARGETS ${library}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
ENDMACRO(SCIRUN_INSTALL_LIBRARY)

# These are installs from the source tree
INSTALL(DIRECTORY pixmaps DESTINATION src
  PATTERN ".svn" EXCLUDE # don't want the subversion directory
  )

INSTALL(FILES 
    scripts/module-remapping.txt 
    DESTINATION src/scripts
  )

INSTALL(FILES
  Applications/SCIRun/scisplash.ppm
  Applications/SCIRun/sciabout.ppm
  Applications/SCIRun/scirun.rc
  DESTINATION src/Applications/SCIRun
  )

INSTALL(FILES
  LICENSE
  DESTINATION bin
  )

INSTALL(FILES
  LICENSE
  DESTINATION src
  )

INSTALL(FILES
  scirunrc
  DESTINATION src
  )

# Installs the directories in src.  Input is a list.  Don't pass stuff
# in as "dir1 dir2 dir3", thi is a single string, not a list.  Pass it
# in as SCIRUN_INSTALL_LIBRARY(dir1 dir2 dir3).
MACRO(SCIRUN_INSTALL_DIRECTORIES)
  FOREACH(dir ${ARGN})
    INSTALL(DIRECTORY ${SCIRUN_CORE_SOURCE_DIR}/${dir}/ DESTINATION src/${dir}
      PATTERN ".svn" EXCLUDE # don't want the subversion directory
      PATTERN "${dir}/*"
      )
  ENDFOREACH(dir)
ENDMACRO(SCIRUN_INSTALL_DIRECTORIES)

# package_name is a value parameter
MACRO(SCIRUN_INSTALL_PACKAGE package_name)
  SET(COPY_DIRS
    Packages/${package_name}/Dataflow/XML
    Packages/${package_name}/Dataflow/GUI
    Packages/${package_name}/nets
    )
  SCIRUN_INSTALL_DIRECTORIES(${COPY_DIRS})

  IF (EXISTS ${SCIRUN_CORE_SOURCE_DIR}/Packages/${package_name}/Services)
    SET(COPY_DIRS
      Packages/${package_name}/Services
      )
    SCIRUN_INSTALL_DIRECTORIES(${COPY_DIRS})
  ENDIF (EXISTS ${SCIRUN_CORE_SOURCE_DIR}/Packages/${package_name}/Services)
ENDMACRO(SCIRUN_INSTALL_PACKAGE)

# Install the core directories
SET(COPY_DIRS
  Dataflow/XML
  Dataflow/GUI
  nets
  )
SCIRUN_INSTALL_DIRECTORIES(${COPY_DIRS})

#################################################################################
# Does several operations for a single library

MACRO(SCIRUN_ADD_LIBRARY library)
  ADD_LIBRARY(${library} ${ARGN})
  IF (WIN32 AND MSVC)
    # Turn off Visual Studio warnings about strcpy et.al. being unsafe
    # for release builds (only supporting CMake default targets).
    SET_TARGET_PROPERTIES(${library}
                         PROPERTIES
                           COMPILE_DEFINITIONS_RELEASE "_CRT_SECURE_NO_WARNINGS=1"
                           COMPILE_DEFINITIONS_MINSIZEREL "_CRT_SECURE_NO_WARNINGS=1"
                         )
  ENDIF (WIN32 AND MSVC)

IF (NOT DISABLE_OPENGL)
  SET_TARGET_PROPERTIES(${library}
                       PROPERTIES
                         COMPILE_DEFINITIONS_DEBUG "OGL_DBG=1")
ENDIF (NOT DISABLE_OPENGL)

  SCIRUN_INSTALL_LIBRARY(${library})

  # Set a temporary variable in order to make sure it stays as a list.
  #SET(temp_library_list ${SCIRUN_LIBRARY_LIST_BUILD} ${library})
  #SET(SCIRUN_LIBRARY_LIST_BUILD ${temp_library_list} CACHE INTERNAL "")
ENDMACRO(SCIRUN_ADD_LIBRARY)

MACRO(SCIRUN_ADD_UNIT_TEST test)
  ADD_EXECUTABLE(${test} ${ARGN})
  LIST(APPEND tmp ${test} ${UNIT_TEST_LIST})
  SET(UNIT_TEST_LIST ${tmp} CACHE INTERNAL "List of unit tests." FORCE)
ENDMACRO(SCIRUN_ADD_UNIT_TEST)

########################################################################
# Windows configuration options

IF(WIN32)
  # FIND_LIBRARY under Windows searches for *.lib files, which are
  # found in the Microsoft SDK.
  #
  # The DLLs for these libraries are found in %SystemRoot%/system32.
  #
  # This registry key should be consistent throughout modern SDK distributions),
  # and we're assuming that we're not going to be building on Itanium.


  IF(NBITS EQUAL 64)
    SET(MS_SDK_REG_PATH
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Microsoft SDKs\\Windows;CurrentInstallFolder]/Lib/x64"
        CACHE INTERNAL
          "Set up search path to Windows SDK."
     )
  ELSE(NBITS EQUAL 64)
    SET(MS_SDK_REG_PATH
          "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Microsoft SDKs\\Windows;CurrentInstallFolder]/Lib"
        CACHE INTERNAL
          "Set up search path to Windows SDK."
     )
  ENDIF(NBITS EQUAL 64)

  IF(MSVC)
    # Bind application to the current version of the VC libraries
    ADD_DEFINITIONS(-D_BIND_TO_CURRENT_VCLIBS_VERSION) 

    # MSVC level 3 warnings are a little ridiculous
    SET(CMAKE_CXX_WARNING_LEVEL 1)
    IF(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
      STRING(REGEX REPLACE "/W[0-4]" "/W1"
              CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    ELSE(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W1")
    ENDIF(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
  ENDIF(MSVC)
	
  #  This defines M_PI for windows
  ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
  # defines Windows 2000 features like mouse leaving window
  ADD_DEFINITIONS(-D_WIN32_WINNT=0x0500)
  ADD_DEFINITIONS(-DWIN32_LEAN_AND_MEAN) # makes sure things like rad1/grp1 don't get defined
  ADD_DEFINITIONS(-DNOMINMAX) # makes sure min and max don't get defined
  
  
  IF(NOT CYGWIN)
    SET(VC "${CMAKE_MAKE_PROGRAM}")
    IF(CMAKE_GENERATOR MATCHES "Visual Studio 7")
      STRING(REGEX REPLACE "Common7/IDE/VCExpress.exe" "Vc7" VC "${VC}")
      STRING(REGEX REPLACE "Common7/IDE/devenv.com" "Vc7" VC "${VC}")
    ELSE(CMAKE_GENERATOR MATCHES "Visual Studio 7")
      STRING(REGEX REPLACE "Common7/IDE/devenv.com" "VC" VC "${VC}")
      STRING(REGEX REPLACE "Common7/IDE/VCExpress.exe" "VC" VC "${VC}")
    ENDIF(CMAKE_GENERATOR MATCHES "Visual Studio 7")
    
    SET (PATH_TO_VC "${VC}" CACHE PATH "Path to Visual Studio")
    MARK_AS_ADVANCED(PATH_TO_VC)

  ENDIF(NOT CYGWIN)
ENDIF(WIN32)

################################################################################################
# SCIRun Thirdparty Externals
#

ADD_DEFINITIONS(-DVPX_CODEC_DISABLE_COMPAT=1)
SET(SCIRUN_EXTERNALS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Externals)
ADD_SUBDIRECTORY(Externals)

IF(NOT BUILD_SHARED_LIBS)
  MESSAGE(STATUS "Building SCIRun Statically")
  ADD_DEFINITIONS(-DBUILD_SCIRUN_STATIC)
  ADD_DEFINITIONS(-DBUILD_THIRDPARTY_STATIC)
ENDIF(NOT BUILD_SHARED_LIBS)

IF(BUILD_STATIC_THIRDPARTY)
  MESSAGE(STATUS "Building Thirdparty Externals Statically")
  ADD_DEFINITIONS(-DBUILD_THIRDPARTY_STATIC)
ENDIF(BUILD_STATIC_THIRDPARTY)

IF(NOT BUILD_SHARED_LIBS OR BUILD_STATIC_THIRDPARTY)
  ADD_DEFINITIONS(-DSLIVR_STATIC)
  ADD_DEFINITIONS(-DGLEW_STATIC)
  # builds against static teem
  ADD_DEFINITIONS(-DTEEM_STATIC)
ELSEIF(NOT BUILD_SHARED_LIBS OR BUILD_STATIC_THIRDPARTY)
  IF(UNIX)
    ADD_DEFINITIONS(-DSLIVR_STATIC)
    ADD_DEFINITIONS(-DGLEW_STATIC)
  ENDIF(UNIX)
ENDIF(NOT BUILD_SHARED_LIBS OR BUILD_STATIC_THIRDPARTY)

########################################################################
# SCIRun Packages to enable (must enable packages when building core for now)
# DV Fix need error checking to make sure SCIRun gets included

SET(LOAD_PACKAGE "SCIRun,BioPSE,Teem,MatlabInterface" CACHE STRING "Building SCIRun, BioPSE, Teem, and MatlabInterface packages by default.")
STRING(REGEX MATCH "MatlabInterface" MATLABINTERFACE_INCLUDED ${LOAD_PACKAGE})
IF (MATLABINTERFACE_INCLUDED)
  ADD_DEFINITIONS(-DWITH_MATLAB)
ENDIF (MATLABINTERFACE_INCLUDED)


########################################################################
# Documentation system

FIND_PACKAGE(LATEX)

OPTION(BUILD_DOCUMENTATION "Build documentation" OFF)

IF (PDFLATEX_COMPILER)
  IF (BUILD_DOCUMENTATION)
    ADD_SUBDIRECTORY(Documentation)
  ENDIF (BUILD_DOCUMENTATION)
ENDIF (PDFLATEX_COMPILER)

########################################################################
# Thread library

FIND_PACKAGE( Threads )
SET(THREAD_LIBRARY ${CMAKE_THREAD_LIBS} CACHE INTERNAL "Thread library.")
# DV Define DEF_THREAD?

IF(CMAKE_USE_PTHREADS)
  FIND_LIBRARY(THREAD_LIBRARY_TEMP pthread)
  IF(THREAD_LIBRARY_TEMP)
    SET(THREAD_LIBRARY ${THREAD_LIBRARY_TEMP} CACHE INTERNAL "Thread library.")
  ENDIF(THREAD_LIBRARY_TEMP)
  
  #SET(DEF_THREAD "#define SCI_PTHREAD 1")
  # Check for sempahore libs w/in pthread library
  # DV can I scale down the lib like configure.ac for loop
  #TRY_COMPILE(HAS_SEMAPHORE_LIBS
  #  ${CMAKE_CURRENT_BINARY_DIR}/CMakeTmp
  #  ${SCIRUN_CORE_SOURCE_DIR}/CMake/testSemaphore.cc
  #    COMPILE_DEFINITIONS
  #    -l${THREAD_LIBRARY}
  #    -I/usr/local/include -I/usr/include -I/include
  #  )
  #IF(HAS_SEMAPHORE_LIBS)
  #  SET(SEMAPHORE_LIBRARY ${THREAD_LIBRARY} CACHE INTERNAL "Semaphore library in pthreads.")
  #ENDIF(HAS_SEMAPHORE_LIBS)

  
ELSE(CMAKE_USE_PTHREADS)
    SET(DEF_THREAD "") 
ENDIF(CMAKE_USE_PTHREADS)

########################################################################
# Socket library
# DV Not sure if socket library needs to be defined for cygwin

IF(WIN32 AND NOT CYGWIN)
  SET(SOCKET_LIBRARY ws2_32 CACHE INTERNAL "Socket library ws2_32")
  SET(THREAD_LIBRARY winmm CACHE INTERNAL "Time functionality library winmm")
ENDIF(WIN32 AND NOT CYGWIN)

########################################################################
# math library

IF(UNIX)
  FIND_LIBRARY(M_LIBRARY m
             PATHS /usr/local/lib
             /usr/lib 
             /lib)
ENDIF(UNIX)

########################################################################
# dl library

IF(UNIX)
  FIND_LIBRARY(DL_LIBRARY dl)
ENDIF(UNIX)

########################################################################
# Search for SGI's or Win32's traceback library
# traceback_lib exc_lib_dir_flag exc_lib_flag
# DV Make required and look for it, also fix for SGI

INCLUDE (${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)
IF(WIN32)
  CHECK_INCLUDE_FILE("windows.h"       HAVE_WINDOWS_H)
  #CHECK_INCLUDE_FILE("dbghlp"       HAVE_DBGHLP_H)
  IF(HAVE_WINDOWS_H)
  #    IF(HAVE_DBGHLP_H)
        SET(TRACEBACK_LIB imagehlp psapi CACHE INTERNAL "Trace back library.")
        SET(HAVE_TRACEBACK 1)
        SET(DEF_TRACEBACK "#define HAVE_TRACEBACK 1")
  #    ELSE(HAVE_DBGHELP_H)
  #      MESSAGE(SEND_ERROR "Trace back library requires dbghlp.h for Windows.")
  #    ENDIF(HAVE_DBGHLP_H)
  ELSE(HAVE_WINDOWS_H)
     MESSAGE(SEND_ERROR "Trace back library requires windows.h for Windows.")
  ENDIF(HAVE_WINDOWS_H)
ELSE(WIN32)
  CHECK_INCLUDE_FILE("libexc.h"       HAVE_TRACEBACK_H)
  IF(HAVE_TRACEBACK_H)
    SET(TRACEBACK_LIB exc CACHE INTERNAL "Trace back library.")
    SET(HAVE_TRACEBACK 1)
    SET(DEF_TRACEBACK "#define HAVE_TRACEBACK 1")
  ENDIF(HAVE_TRACEBACK_H)
ENDIF(WIN32)

########################################################################
# OpenGL library - win32 finds the OPENGL_gl_LIBRARY but for some reason this 
#   if check does not complete
# Use our own find OpenGL module

SET (CMAKE_MODULE_PATH ${SCIRUN_CORE_SOURCE_DIR}/CMake)
FIND_PACKAGE( OpenGL )
IF(OPENGL_LIBRARIES AND OPENGL_INCLUDE_DIR)
  SET(GL_LIBRARY ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} CACHE STRING "OpenGL Libs")
  SET(GL_INCLUDE ${OPENGL_INCLUDE_DIR} CACHE STRING "OpenGL Includes")
  SET(SCIRUN_INCLUDE_ADDITIONAL ${SCIRUN_INCLUDE_ADDITIONAL} ${GL_INCLUDE})
ENDIF(OPENGL_LIBRARIES AND OPENGL_INCLUDE_DIR)

IF (WIN32)
  SET(GL_LIBRARY ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} gdi32 CACHE STRING "OpenGL Libs")
  SET(GL_INCLUDE ${OPENGL_INCLUDE_DIR} CACHE STRING "OpenGL Includes")
  SET(SCIRUN_INCLUDE_ADDITIONAL ${SCIRUN_INCLUDE_ADDITIONAL} ${GL_INCLUDE})
ENDIF(WIN32)


########################################################################
# X11
OPTION(WITH_X11 "Build SCIRun with X11 support" ON)

SET(X11_FOUND 0)

IF(WITH_X11)
  FIND_PACKAGE( X11 )
  IF(X11_FOUND)
    SET(X_LIBRARY ${X11_LIBRARIES})
    SET(DEF_X11 "#define HAVE_X11 1")
  ENDIF(X11_FOUND)

  IF(APPLE)
  ADD_DEFINITIONS(-DGLEW_APPLE_GLX)
  ENDIF(APPLE)
ENDIF(WITH_X11)

MARK_AS_ADVANCED(WITH_X11)

# Fixing a problem on OSX for now: the standard CMake script finds the wrong OpenGL
# libraries. Hence they are hard coded here, however these require X11 to be installed
# and this is not checked here. We need to come up with a better fix here.

IF (APPLE AND X11_FOUND)

  IF (APPLE_LEOPARD)
    SET (GL_LIBRARY "" CACHE INTERNAL "OpenGL Libs" FORCE)
    SET (X_LIBRARY "" CACHE INTERNAL "OpenGL Libs" FORCE)
    IF (NBITS EQUAL 32)
      SET (GL_LINK "-dylib_file /System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib /usr/X11R6/lib/libGL.dylib /usr/X11R6/lib/libGLU.dylib  -framework AGL" CACHE INTERNAL "Link options for OpenGL libs" FORCE)
      SET (X_LINK " /usr/X11R6/lib/libX11.dylib /usr/X11R6/lib/libXext.dylib -framework Carbon" CACHE INTERNAL "Link options for X11 libs" FORCE)
    ELSE (NBITS EQUAL 32)
      SET (GL_LINK "-dylib_file /System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib /usr/X11R6/lib/libGL.dylib /usr/X11R6/lib/libGLU.dylib  -framework AGL" CACHE INTERNAL "Link options for OpenGL libs" FORCE)
      SET (X_LINK " /usr/X11R6/lib/libX11.dylib /usr/X11R6/lib/libXext.dylib" CACHE INTERNAL "Link options for X11 libs" FORCE)  
    ENDIF (NBITS EQUAL 32)
  
  ELSE (APPLE_LEOPARD)
    SET (GL_LIBRARY /usr/X11R6/lib/libGL.dylib /usr/X11R6/lib/libGLU.dylib "-framework AGL" CACHE INTERNAL "OpenGL Libs")
    SET (GL_LINK "" CACHE INTERNAL "Link options for OpenGL libs" FORCE)
    SET (X_LINK "" CACHE INTERNAL "Link options for X libs" FORCE)
    SET (SCIRUN_INCLUDE_ADDITIONAL ${SCIRUN_INCLUDE_ADDITIONAL} ${GL_INCLUDE})    
  ENDIF (APPLE_LEOPARD)

ENDIF(APPLE AND X11_FOUND)

########################################################################
# sci_defs variables

# sci_defs/compile_defs.h
# DV Not needed eventually?

# Set Make for MAKE_COMMAND in compile_defs.h
# DV Fix Make sure this works for win
SET(MAKE ${CMAKE_BUILD_TOOL})

# sci_defs/environment_defs.h
SET(scirun_src_dir "${SCIRUN_CORE_SOURCE_DIR}" CACHE INTERNAL "Path to SCIRun's src directory")
# DV Technically, I don't have access to this variable once we fix building
# these projets separately!!!
SET(scirun_obj_dir ${SCIRUN_CORE_BINARY_DIR} CACHE INTERNAL "Path to SCIRun's binary directory")

# sci_defs/error_defs.h
SET(THROW_LEVEL 0 CACHE STRING "Set throw level N(0-3)")
SET(ASSERTION_LEVEL 1 CACHE STRING "Set assertion level N(0-3)")
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  SET(ASSERTION_LEVEL 3)
ENDIF(CMAKE_BUILD_TYPE MATCHES "Debug")


# sci_defs/lapack_defs.h
# DV FIX hard coding of lapack library
IF(HAVE_LAPACK)
  SET(DEF_LAPACK "#define HAVE_LAPACK 1" CACHE INTERNAL "Have LAPACK")
ELSE(HAVE_LAPACK)
  UNSET(DEF_LAPACK CACHE)
ENDIF(HAVE_LAPACK)

IF(HAVE_CBLAS)
  SET(DEF_CBLAS "#define HAVE_CBLAS 1" CACHE INTERNAL "Have CBLAS")
ELSE(HAVE_CBLAS)
  UNSET(DEF_CBLAS CACHE)
ENDIF(HAVE_CBLAS)

IF(HAVE_BLAS)
  SET(DEF_BLAS "#define HAVE_BLAS 1" CACHE INTERNAL "Have BLAS")
ELSE(HAVE_BLAS)
  UNSET(DEF_BLAS CACHE)
ENDIF(HAVE_BLAS)

IF(BUILD_NETLIB)
  SET(DEF_SCIRUN_NETLIB "#define HAVE_SCIRUN_NETLIB 1" CACHE INTERNAL "Use internal netlib LAPACK library")
ELSE(BUILD_NETLIB)
  UNSET(DEF_SCIRUN_NETLIB CACHE)
ENDIF(BUILD_NETLIB)

# TODO: replace with Boost unordered map
#
# sci_defs/hashmap_defs.h
SET(DEF_HASHMAP "None" CACHE INTERNAL "Determining hash map")

TRY_COMPILE(DEF_HASHMAP
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeTmp
  ${SCIRUN_CORE_SOURCE_DIR}/CMake/SCIRunTestUnorderedSet.cc
  OUTPUT_VARIABLE OUTPUT)
IF(DEF_HASHMAP)
  SET(DEF_HASHMAP "#define HAVE_STD_UNORDERED_SET 1" CACHE INTERNAL
      "Using STD unordered set.")
ELSE(DEF_HASHMAP)
  # DV Fix to figure out value of DEF_HASHMAP
  TRY_COMPILE(DEF_HASHMAP
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeTmp
    ${SCIRUN_CORE_SOURCE_DIR}/CMake/SCIRunTestHashMapSTD.cc
    OUTPUT_VARIABLE OUTPUT)
  IF(DEF_HASHMAP)
    SET(DEF_HASHMAP "#define HAVE_STD_HASHMAP 1" CACHE INTERNAL "Using STD hash map")
  ELSE(DEF_HASHMAP)
    TRY_COMPILE(DEF_HASHMAP
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeTmp
    ${SCIRUN_CORE_SOURCE_DIR}/CMake/SCIRunTestHashMapEXT.cc
    OUTPUT_VARIABLE OUTPUT)

    # include <ext/hash_map> instead
    IF(DEF_HASHMAP)
      SET(DEF_HASHMAP "#define HAVE_EXT_HASHMAP 1" CACHE INTERNAL "Using EXT hash map")
    ELSE(DEF_HASHMAP)
     TRY_COMPILE(DEF_HASHMAP
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeTmp
      ${SCIRUN_CORE_SOURCE_DIR}/CMake/SCIRunTestHashMapSTDEXT.cc
      OUTPUT_VARIABLE OUTPUT)

      IF(DEF_HASHMAP)
        SET(DEF_HASHMAP "#define HAVE_STDEXT_HASHMAP 1" CACHE INTERNAL
            "Using STD EXT hash map")
      ELSE(DEF_HASHMAP)
        TRY_COMPILE(DEF_HASHMAP
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeTmp
        ${SCIRUN_CORE_SOURCE_DIR}/CMake/SCIRunTestHashMapGNU.cc
        OUTPUT_VARIABLE OUTPUT)
        # include <ext/hash_map> instead

        IF(DEF_HASHMAP)
          SET(DEF_HASHMAP "#define HAVE_GNU_HASHMAP 1" CACHE INTERNAL "Using GNU hash map")
        ELSE(DEF_HASHMAP)
          SET(DEF_HASHMAP "" CACHE INTERNAL "No hash map found")
          MESSAGE(STATUS "Cannot find hash_map, will use map instead")
        ENDIF(DEF_HASHMAP)
      ENDIF(DEF_HASHMAP)
    ENDIF(DEF_HASHMAP)
  ENDIF(DEF_HASHMAP)
ENDIF(DEF_HASHMAP)

MARK_AS_ADVANCED(THROW_LEVEL ASSERTION_LEVEL DEF_HASH)
# DV Add with_exc stuff from configure.ac
INCLUDE (${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)
CHECK_INCLUDE_FILE(limits.h HAVE_LIMITS)
IF(HAVE_LIMITS)
  SET(HAVE_LIMITS_H 1)
ELSE(HAVE_LIMITS)
  SET(HAVE_LIMITS_H 0)
ENDIF(HAVE_LIMITS)

CHECK_INCLUDE_FILE(sys/select.h HAVE_SYS_SELECT)
IF(HAVE_SYS_SELECT)
  SET(HAVE_SYS_SELECT_H 1)
ELSE(HAVE_SYS_SELECT)
  SET(HAVE_SYS_SELECT_H 0)
ENDIF(HAVE_SYS_SELECT)

CHECK_INCLUDE_FILE(sys/time.h HAVE_SYS_TIME)
IF(HAVE_SYS_TIME)
  SET(HAVE_SYS_TIME_H 1)
ELSE(HAVE_SYS_TIME)
  SET(HAVE_SYS_TIME_H 0)
ENDIF(HAVE_SYS_TIME)

# BW Fix to figure out if we have ext/algorithm
TRY_COMPILE(HAVE_EXT_ALGORITHM
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeTmp
  ${SCIRUN_CORE_SOURCE_DIR}/CMake/SCIRunTestExtAlgorithm.cc
  OUTPUT_VARIABLE OUTPUT)
IF(HAVE_EXT_ALGORITHM)
  SET(HAVE_EXT_ALGORITHM_H 1)
ELSE(HAVE_EXT_ALGORITHM)
  SET(HAVE_EXT_ALGORITHM_H 0)
ENDIF(HAVE_EXT_ALGORITHM)

CHECK_INCLUDE_FILE(stdint.h HAVE_STDINT)
IF(HAVE_STDINT)
  SET(HAVE_STDINT_H 1)
ELSE(HAVE_STDINT)
  SET(HAVE_STDINT_H 0)
ENDIF(HAVE_STDINT)

CHECK_INCLUDE_FILE(inttypes.h HAVE_INTTYPES)
IF(HAVE_INTTYPES)
  SET(HAVE_INTTYPES_H 1)
ELSE(HAVE_INTTYPES)
  SET(HAVE_INTTYPES_H 0)
ENDIF(HAVE_INTTYPES)

SET(HAVE_SYS_UUID_H 0)
SET(HAVE_UUID_UUID_H 0)


########################################################################
# hypre library
# sci_defs/hypre_defs.h

OPTION(WITH_HYPRE "Build with hypre library" OFF)
MARK_AS_ADVANCED(WITH_HYPRE)		   
IF(WITH_HYPRE OR HYPRE_DIR)
  MESSAGE(STATUS "Looking for Hypre")
  SET(HYPRE_DIR ${HYPRE_DIR} CACHE PATH "Path to Hypre" FORCE)
  MARK_AS_ADVANCED(CLEAR HYPRE_DIR)

  FIND_LIBRARY(HYPRE_LIBRARY NAMES HYPRE PATHS ${HYPRE_DIR}/lib)
  FIND_PATH(HYPRE_INCLUDE HYPRE.h PATHS ${HYPRE_DIR}/include)

  IF(HYPRE_LIBRARY AND HYPRE_INCLUDE)
    MARK_AS_ADVANCED(HYPRE_LIBRARY)
    MARK_AS_ADVANCED(HYPRE_INCLUDE)
    MESSAGE(STATUS "Looking for hypre - found")
    SET(HAVE_HYPRE "YES")
    SET(DEF_HYPRE "#define HAVE_HYPRE 1")
  ELSE(HYPRE_LIBRARY AND HYPRE_INCLUDE)
    MESSAGE(STATUS "Looking for hypre - not found")
    SET(HAVE_HYPRE "NO")
  ENDIF(HYPRE_LIBRARY AND HYPRE_INCLUDE)
ELSE(WITH_HYPRE OR HYPRE_DIR)
   SET(HAVE_HYPRE "NO")
ENDIF(WITH_HYPRE OR HYPRE_DIR)

# sci_defs/version_defs.h

IF(WIN32)
  SET(CMAKE_PROGRAM_PATH_OLD ${CMAKE_PROGRAM_PATH})
  SET(CMAKE_PROGRAM_PATH "c:/cygwin/bin")
ENDIF(WIN32)

FIND_FILE(REVISION_FILE REVISION.txt 
          PATHS ${SCIRUN_CORE_SOURCE_DIR}
          NO_DEFAULT_PATH)

IF(REVISION_FILE)
  FILE(READ ${REVISION_FILE} SCIRUN_SVN_REVISION)
  SET(SCIRUN_SVN_URL "https://code.sci.utah.edu/svn/SCIRun/cibc/trunk/SCIRun/src")
ELSE(REVISION_FILE)
FIND_PACKAGE(Subversion)
IF(Subversion_FOUND)
  Subversion_WC_INFO(${SCIRUN_CORE_SOURCE_DIR} SCIRUN_SVN)
ENDIF(Subversion_FOUND)
ENDIF(REVISION_FILE)

IF(WIN32)
  SET(CMAKE_PROGRAM_PATH ${CMAKE_PROGRAM_PATH_OLD})
ENDIF(WIN32)

########################################################################
# Insight Toolkit

STRING(REGEX MATCH "Insight" BUILDING_INSIGHT ${LOAD_PACKAGE})
IF(WITH_ITK)
  # Check if Insight in LOAD_PACKAGE
  IF(NOT BUILDING_INSIGHT)
    MESSAGE(STATUS "Building ITK, but Insight not included in LOAD_PACKAGE. Please add Insight to LOAD_PACKAGE if package functionality is needed.")
    MESSAGE(STATUS "LoadPackage=${LOAD_PACKAGE}")
  ENDIF(NOT BUILDING_INSIGHT)
ELSE(WITH_ITK)
  # Force user to define it if they specified Insight in LOAD_PACKAGE
  IF(BUILDING_INSIGHT)
    MESSAGE(SEND_ERROR "Insight specified in LOAD_PACKAGE, but WITH_ITK not set. Please set WITH_ITK to ON or remove Insight from LOAD_PACKAGE")
  ENDIF(BUILDING_INSIGHT)
ENDIF(WITH_ITK)

IF(WITH_ITK)
    SET(DEF_INSIGHT "#define HAVE_INSIGHT 1" CACHE INTERNAL "Have insight")
    SET(HAVE_INSIGHT "yes" CACHE INTERNAL "Linking against ITK")
  SET(ITK_FOUND 1 CACHE INTERNAL "ITK has been found")	
ENDIF(WITH_ITK)

########################################################################
# Search for HDF5 (http://www.hdfgroup.org)

# TODO: write a module for this

IF(HDF5_DIR)
  SET(HDF5_DIR ${HDF5_DIR} CACHE PATH "Path to HDF5" FORCE)
  MARK_AS_ADVANCED(CLEAR HDF5_DIR)
ENDIF(HDF5_DIR)

IF(HDF5_DIR)

  MESSAGE(STATUS "Configuring SCIRun with HDF5")

  SET(HDF5_LIBRARY HDF5_LIBRARY-NOTFOUND CACHE FILEPATH 
      "HDF5 library" FORCE)
  SET(HDF5_INCLUDE HDF5_INCLUDE-NOTFOUND CACHE FILEPATH 
      "HDF5 include" FORCE)

	IF(HDF5_DIR)
    FIND_LIBRARY(HDF5_LIBRARY hdf5 
			${HDF5_DIR}
			${HDF5_DIR}/lib)
    FIND_PATH(HDF5_INCLUDE hdf5.h 
			${HDF5_DIR}
			${HDF5_DIR}/include
			${HDF5_DIR}/../include)
	ENDIF(HDF5_DIR)
	
  IF(HDF5_LIBRARY AND HDF5_INCLUDE)
    SET(HAVE_HDF5 "YES" CACHE INTERNAL "Do not edit" FORCE)
    SET(DEF_HDF5 "#define HAVE_HDF5 1" CACHE INTERNAL "Do not edit" FORCE)
    SET(SCIRUN_INCLUDE_THIRDPARTY ${SCIRUN_INCLUDE_THIRDPARTY} 
        ${HDF5_INCLUDE})
  ELSE(HDF5_LIBRARY AND HDF5_INCLUDE)
    MESSAGE(FATAL_ERROR "Cannot find HDF5 library.") 
  ENDIF(HDF5_LIBRARY AND HDF5_INCLUDE)

  MARK_AS_ADVANCED(HDF5_LIBRARY)
  MARK_AS_ADVANCED(HDF5_INCLUDE)
ELSE(HDF5_DIR)
  SET(HAVE_HDF5 "NO" CACHE INTERNAL "Do not edit" FORCE)
  SET(DEF_HDF5 "" CACHE INTERNAL "Do not edit" FORCE)
ENDIF(HDF5_DIR)

########################################################################
# GDCM Library (which could be pulled from ITK or specified individually)
IF(WITH_GDCM)
  # If building ITK, all of the following will be set appropriately
  # by Externals/itk/CMakeLists.txt based on ITK_USE_SYSTEM_GDCM
  # Only set these here if not building ITK
  IF(NOT WITH_ITK)
    SET(GDCM_INCLUDE_DIR "${GDCM_DIR}/src" CACHE PATH "Location of gdcm.h")
    SET(GDCM_CONFIGURE_INCLUDE_DIR "${GDCM_DIR}/bin" CACHE PATH "Location of gdcmConfigure.h")
    SET(GDCM_LIBRARY_DIR "${GDCM_DIR}/bin" CACHE PATH "Location of libgdcm or libitkgdcm")
    MARK_AS_ADVANCED(GDCM_INCLUDE_DIR GDCM_CONFIGURE_INCLUDE_DIR GDCM_LIBRARY_DIR)
    FIND_LIBRARY(GDCM_LIBRARY NAMES gdcm itkgdcm
      PATHS ${GDCM_LIBRARY_DIR} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/InsightToolkit
    ${ITK_DIR}/bin
    PATHS ${SCIRUN_CORE_BINARY_DIR}/lib/gdcm
    NO_DEFAULT_PATH)
  FIND_PATH(GDCM_INCLUDE gdcm.h
      PATHS ${GDCM_INCLUDE_DIR} ${SCIRUN_CORE_BINARY_DIR}/include/InsightToolkit/gdcm/src/
	${ITK_DIR}/include/InsightToolkit/gdcm/src/
    ${ITK_DIR}/Utilities/gdcm/src
    ${SCIRUN_CORE_BINARY_DIR}/include/gdcm
    NO_DEFAULT_PATH)
  FIND_PATH(GDCM_CONFIGURE_INCLUDE gdcmConfigure.h
      PATHS ${GDCM_CONFIGURE_INCLUDE_DIR} ${SCIRUN_CORE_BINARY_DIR}/include/InsightToolkit/gdcm/src/
	${ITK_DIR}/include/InsightToolkit/gdcm/src
    ${ITK_DIR}/Utilities/gdcm
    ${SCIRUN_CORE_BINARY_DIR}/include/gdcm
    NO_DEFAULT_PATH)
  ENDIF(NOT WITH_ITK)
  IF(GDCM_LIBRARY)
    SET(SCIRUN_INCLUDE_THIRDPARTY ${SCIRUN_INCLUDE_THIRDPARTY} ${GDCM_INCLUDE} ${GDCM_CONFIGURE_INCLUDE})
    SET(DEF_GDCM "#define HAVE_GDCM 1")
  ELSE(GDCM_LIBRARY)
    SET(DEF_GDCM "#define HAVE_GDCM 0")
  ENDIF(GDCM_LIBRARY)
ENDIF(WITH_GDCM)

MARK_AS_ADVANCED(WITH_GDCM GDCM_INCLUDE GDCM_CONFIGURE_INCLUDE GDCM_LIBRARY)


#############################################################################################
# The entire SCIRun tree should use the same include path
#
# TODO: amend these comments
# This should be the only INCLUDE_DIRECTORIES command in the entire
# tree.  We need to do this in one place to make sure the order is correct.

# Include directories for other projects installed on the system.
SET(SCIRUN_INCLUDE_DIRS_SYSTEM 
  ${SCIRUN_CORE_SOURCE_DIR}/
  ${SCIRUN_CORE_SOURCE_DIR}/include/
)

# Set include directories for sci_defs

## This has been modified temporarily for new Externals config
## and should be revised to be more modular

SET(SCIRUN_INCLUDE_DIRS_BUILD_TREE 
  ${CMAKE_CURRENT_BINARY_DIR}/include/
  ${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs/
)

MESSAGE(STATUS "${SCI_INCLUDE_EXTERNAL}")

INCLUDE_DIRECTORIES(
  ${SCI_ZLIB_INCLUDE}
  ${SCI_PNG_INCLUDE}
  ${SCI_TEEM_INCLUDE}
  ${SCI_CLAPACK_INCLUDE}
  ${SCI_GLEW_INCLUDE}
  ${SCI_SLIVR_INCLUDE}
  ${SCI_FREETYPE_INCLUDE}
  ${SCI_ITK_INCLUDE} 
  ${SCI_TCL_INCLUDE} 
  ${SCI_TK_INCLUDE} 
  ${SCI_ITCL_INCLUDE} 
  ${SCI_IncrTK_INCLUDE} 
  ${SCI_LIBXML2_INCLUDE}
  ${SCI_TETGEN_INCLUDE}
  ${SCI_BOOST_INCLUDE}
  ${SCI_GTEST_INCLUDE}
  ${SCIRUN_INCLUDE_DIRS_BUILD_TREE}
  ${SCIRUN_INCLUDE_DIRS_SYSTEM}
  ${SCIRUN_INCLUDE_ADDITIONAL}
  ${SCIRUN_INCLUDE_THIRDPARTY}
  ${SCI_VPX_INCLUDE}
)

########################################################################
# for now set make install directory to current binary directory

SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "Install path")

# StandAlone apps
OPTION(BUILD_UTILS "Build Command line utilities" ON)

IF(BUILD_UTILS)
  ADD_SUBDIRECTORY(Applications/Utils)
ENDIF(BUILD_UTILS)


################################################################################################
# Option for building scirun executable, dataflow, packages
IF(BUILD_DATAFLOW)
  SET(SCIRUN_DIR ${SCIRUN_CORE_BINARY_DIR})
  MARK_AS_ADVANCED(SCIRUN_DIR)
  ADD_SUBDIRECTORY(Dataflow)
  ADD_SUBDIRECTORY(Applications/SCIRun)
ENDIF(BUILD_DATAFLOW)

# Traverse Testing directory
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(Testing)
  SET(UNIT_TEST_LIST "" CACHE INTERNAL "List of unit tests." FORCE)
  ADD_DEFINITIONS(-DUNIT_TEST_DATA_ROOT="${SCIRUN_UNIT_TEST_DATA_DIR}")
ENDIF(BUILD_TESTING)

MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX DART_TESTING_TIMEOUT)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)

# traverse the tree
ADD_SUBDIRECTORY(include)
ADD_SUBDIRECTORY(Core)

IF (NOT WIN32)
  ADD_CUSTOM_TARGET(copy_scripts ALL
    COMMAND cp ${SCIRUN_CORE_SOURCE_DIR}/scripts/run_unit_tests.sh ${SCIRUN_CORE_BINARY_DIR}/run_unit_tests.sh
    COMMAND chmod +x ${SCIRUN_CORE_BINARY_DIR}/run_unit_tests.sh
    COMMAND cp ${SCIRUN_CORE_SOURCE_DIR}/scripts/run_sample_tests.sh ${SCIRUN_CORE_BINARY_DIR}/run_sample_tests.sh
    COMMAND chmod +x ${SCIRUN_CORE_BINARY_DIR}/run_sample_tests.sh
    COMMAND cp ${SCIRUN_CORE_SOURCE_DIR}/scripts/run_all_tests.sh ${SCIRUN_CORE_BINARY_DIR}/run_all_tests.sh
    COMMAND chmod +x ${SCIRUN_CORE_BINARY_DIR}/run_all_tests.sh
#    COMMENT "Copying scripts for testing/using dashboard"
    )
ENDIF(NOT WIN32)

# This needs to happen after all the directories are traversed.  The following
# is a hack to get around an annoying scoping issue.  SCIRUN_LIBRARY_LIST needs
# to start out blank in order to get built fresh during cmake invocation.  If
# you set it blank at the start of the run, then it doesn't build up properly
# (it basically results in getting the last value set instead of concatenation).
# The fix is to add an additional variable, SCIRUN_LIBRARY_LIST_BUILD.  This
# variable is built up each time, is copied to the real variable, then erased
# for the next run.

SET(SCIRUN_LIBRARY_LIST ${SCIRUN_LIBRARY_LIST_BUILD} CACHE INTERNAL "List of all libraries")
SET(SCIRUN_LIBRARY_LIST_BUILD "" CACHE INTERNAL "Configure time list of all libraries (should be blank when done)")

# MESSAGE("SCIRUN_LIBRARY_LIST = ")
# FOREACH(library ${SCIRUN_LIBRARY_LIST})
#   MESSAGE("lib = ${library}")
# ENDFOREACH(library)

# Some stuff for installing stuff
ADD_SUBDIRECTORY(InstallUtils)

# License files

FILE(MAKE_DIRECTORY "${SCIRUN_CORE_BINARY_DIR}/Licenses")

CONFIGURE_FILE(${SCIRUN_CORE_SOURCE_DIR}/LICENSE
                           ${SCIRUN_CORE_BINARY_DIR}/Licenses/SCIRun.license
                           COPYONLY
                           )

IF(SCI_TEEM_LIBRARY)
  CONFIGURE_FILE(${SCIRUN_CORE_SOURCE_DIR}/Externals/teem/LICENSE.txt
                             ${SCIRUN_CORE_BINARY_DIR}/Licenses/Teem.license
                              COPYONLY
                            )
ENDIF(SCI_TEEM_LIBRARY)

IF(SCI_TETGEN_LIBRARY)
  CONFIGURE_FILE(${SCIRUN_CORE_SOURCE_DIR}/Externals/tetgen/LICENSE
                             ${SCIRUN_CORE_BINARY_DIR}/Licenses/TetGen.license
                             COPYONLY
                            )
ENDIF(SCI_TETGEN_LIBRARY)


