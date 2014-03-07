#Then for each subdirectory listed in the SUBDIRS command,
# CMakeLists.txt files are created. In the ./Hello directory, the
# following CMakeLists.txt file is created:

# Create a library called "Hello" which includes the source file "hello.cxx".
# Any number of sources could be listed here.
ADD_LIBRARY(util array.cxx mathutil.cxx list.cxx)

IF(WIN32)
  SET_TARGET_PROPERTIES(util
                        PROPERTIES
                          COMPILE_FLAGS "-DBUILD_util"
                          DEFINE_SYMBOL _USRDLL
                       )
ENDIF(WIN32)

