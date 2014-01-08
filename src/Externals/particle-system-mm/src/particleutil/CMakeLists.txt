SET(PARTICLEUTIL_SRCS
       NrrdIO.cxx
      )

SET(PARTICLEUTIL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${TEEM_INCLUDE_DIRS} CACHE INTERNAL "Particle system utility classes.")

INCLUDE_DIRECTORIES(${PARTICLEUTIL_INCLUDES})

ADD_LIBRARY(particleutil ${PARTICLEUTIL_SRCS})

TARGET_LINK_LIBRARIES(particleutil ${TEEM_LIBRARIES})
IF(WIN32)
  SET_TARGET_PROPERTIES(particleutil
                        PROPERTIES
                         COMPILE_FLAGS "-DBUILD_Particleutil"
                         DEFINE_SYMBOL _USRDLL
                       )
ENDIF(WIN32)
