SET(UTIL_SRCS
  IOScalarField.cxx
  mtxlib.cxx
  ScalarField.cxx
  Surface.cxx
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${PARTICLEUTIL_INCLUDES})

ADD_LIBRARY(mm_util ${UTIL_SRCS})
TARGET_LINK_LIBRARIES(mm_util particleutil)
IF(WIN32)
  SET_TARGET_PROPERTIES(mm_util
                        PROPERTIES
                          COMPILE_FLAGS "-DBUILD_MM_util"
                          DEFINE_SYMBOL _USRDLL
                       )
ENDIF(WIN32)

SET(UTIL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${PARTICLEUTIL_INCLUDES})
