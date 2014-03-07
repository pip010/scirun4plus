#Then for each subdirectory listed in the SUBDIRS command,
# CMakeLists.txt files are created. In the ./Hello directory, the
# following CMakeLists.txt file is created:

# Create a library called "Hello" which includes the source file "hello.cxx".
# Any number of sources could be listed here.
ADD_LIBRARY(matrix matrix.cxx vector.cxx)

IF(WITH_NETLIB)
  ADD_DEFINITIONS("-DUSELAPACKLIB")
ENDIF(WITH_NETLIB)

IF(WIN32)
  SET_TARGET_PROPERTIES(matrix
                        PROPERTIES
                          COMPILE_FLAGS "-DBUILD_matrix"
                          DEFINE_SYMBOL _USRDLL
                       )
ENDIF(WIN32)
