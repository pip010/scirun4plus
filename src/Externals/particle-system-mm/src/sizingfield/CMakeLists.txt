SET(LFS_SRCS
  LFS.cxx 
  SizingField.cxx
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${PARTICLEUTIL_INCLUDES})

ADD_EXECUTABLE(sizingfield ${LFS_SRCS} main.cpp)
TARGET_LINK_LIBRARIES(sizingfield mm_util particleutil)

ADD_EXECUTABLE(redosizingfield ${LFS_SRCS} main_redo.cpp)
TARGET_LINK_LIBRARIES(redosizingfield mm_util particleutil)
