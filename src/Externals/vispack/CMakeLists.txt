#The first, top-level directory contains the following CMakeLists.txt file.
# The name of our project is "VISPACK".  CMakeLists files in this project can
# refer to the root source directory of the project as ${VISPACK_SOURCE_DIR} and
# to the root binary directory of the project as ${VISPACK_BINARY_DIR}.

# Vispack build configurations:
#
# BUILD_BIOMESH3D is set in the top-level SCIRun CMakeLists.txt file
# and is used when building applications used in the BioMesh3D pipeline.
#
# BUILD_VISPACK_LAPACK is set in the top-level SCIRun CMakeLists.txt
# file and is used when only building the lapack package.

IF(BUILD_BIOMESH3D OR BUILD_VISPACK_LAPACK)

 CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
  IF(COMMAND CMAKE_POLICY)
    CMAKE_POLICY(VERSION 2.6)
  ENDIF(COMMAND CMAKE_POLICY)
  
  INCLUDE(CMakeBackwardCompatibilityC)

ELSE(BUILD_BIOMESH3D OR BUILD_VISPACK_LAPACK)

  CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

ENDIF(BUILD_BIOMESH3D OR BUILD_VISPACK_LAPACK)

PROJECT(VISPACK)

IF(BUILD_VISPACK_LAPACK)

  SUBDIRS(packages)

ELSE(BUILD_VISPACK_LAPACK)

  IF(UNIX)
    FIND_LIBRARY(M_LIBRARY m)
  ENDIF(UNIX)

  # Recurse into the "Hello" and "Demo" subdirectories.  This does not actually
  # cause another cmake executable to run.  The same process will walk through
  # the project's entire directory structure.
  #SUBDIRS(src packages apps doc)
  SUBDIRS(packages src apps)

  SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${VISPACK_SOURCE_DIR})

  # The BioMesh3D pipeline uses Teem configured for the SCIRun build.
  IF(BUILD_BIOMESH3D)
    SET(TEEM_INCLUDE_DIRS ${SCI_TEEM_INCLUDE})
    SET(TEEM_LIBRARIES ${SCI_TEEM_LIBRARY} ${SCI_PNG_LIBRARY} ${SCI_ZLIB_LIBRARY})
    SET(TEEM_LIBRARY_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

    IF(CMAKE_HAVE_UNISTD_H)
      ADD_DEFINITIONS(-DHAVE_UNISTD_H)
    ENDIF(CMAKE_HAVE_UNISTD_H)

  ELSE(BUILD_BIOMESH3D)
    FIND_PACKAGE(TEEM REQUIRED)
  ENDIF(BUILD_BIOMESH3D)

  SET(VISPACK_PACKAGES_DIR "${VISPACK_SOURCE_DIR}/packages")
  SET(LIBRARY_OUTPUT_PATH ${VISPACK_BINARY_DIR}/lib)

  IF(BUILD_BIOMESH3D)
    SET(VISPACK_INCLUDE_DIRECTORIES
          ${VISPACK_SOURCE_DIR}/src
          ${VISPACK_SOURCE_DIR}/src/util
          ${VISPACK_SOURCE_DIR}/src/image
          ${VISPACK_SOURCE_DIR}/src/vol
          ${VISPACK_SOURCE_DIR}/src/mat
          ${VISPACK_SOURCE_DIR}/src/param
   )

   IF(WIN32 AND MSVC)
     # Bind application to the current version of the VC libraries
     ADD_DEFINITIONS(-D_BIND_TO_CURRENT_VCLIBS_VERSION) 

     # MSVC level 3 warnings are a little ridiculous
     SET(CMAKE_CXX_WARNING_LEVEL 1)
     SET(CMAKE_C_WARNING_LEVEL 1)
     IF(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
       STRING(REGEX REPLACE "/W[0-4]" "/W1"
                CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
     ELSE(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
       SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W1")
     ENDIF(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")

     IF(CMAKE_C_FLAGS MATCHES "/W[0-4]")
       STRING(REGEX REPLACE "/W[0-4]" "/W1"
                CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
     ELSE(CMAKE_C_FLAGS MATCHES "/W[0-4]")
       SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W1")
     ENDIF(CMAKE_C_FLAGS MATCHES "/W[0-4]")
   ENDIF(WIN32 AND MSVC)

  ELSE(BUILD_BIOMESH3D)
    SET(VISPACK_INCLUDE_DIRECTORIES
        ${VISPACK_SOURCE_DIR}/src
        ${VISPACK_SOURCE_DIR}/src/util
        ${VISPACK_SOURCE_DIR}/src/image
        ${VISPACK_SOURCE_DIR}/src/vol
        ${VISPACK_SOURCE_DIR}/src/mat
        ${VISPACK_SOURCE_DIR}/src/voxmodel
        ${VISPACK_SOURCE_DIR}/src/scan
        ${VISPACK_SOURCE_DIR}/src/param
        ${VISPACK_SOURCE_DIR}/src/diffuse
    )
  ENDIF(BUILD_BIOMESH3D)

  #Finally, in the ./Demo directory, the third and final CMakeLists.txt
  #file is created: 
  # Make sure the compiler can find include files from our Hello library.
  INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES} ${TEEM_INCLUDE_DIRS})
  #INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES})

  #Need this for the sgi
  #SET(CMAKE_CXX_FLAGS -LANG:std)
  #SET(CMAKE_CXX_FLAGS -O2)

ENDIF(BUILD_VISPACK_LAPACK)
