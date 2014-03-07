SET(features_SRCS
  mtxlib.cpp
)

ADD_LIBRARY(features ${features_SRCS})
IF(WIN32)
  SET_TARGET_PROPERTIES(features
                        PROPERTIES
                          COMPILE_FLAGS "-DBUILD_features"
                          DEFINE_SYMBOL _USRDLL
                       )
ENDIF(WIN32)
