IF(BUILD_BIOMESH3D)

  # Make sure the compiler can find include files from our Hellolibrary.
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${VISPACK_INCLUDE_DIRECTORIES} ${NRRDIO_INCLUDE})

  # Make sure the linker can find the Hello library once it is built.
  LINK_DIRECTORIES(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${NRRDIO_LIB})

  ADD_EXECUTABLE(linear_test linear_test.cxx)
  TARGET_LINK_LIBRARIES(linear_test
    matrix
    util
    lapack2
    lapack1
    blas
    F77
    ${M_LIBRARY})

ELSE(BUILD_BIOMESH3D)

  #SET(CMAKE_CXX_FLAGS -O3)
  SET(CMAKE_CXX_FLAGS -g)

  # Make sure the compiler can find include files from our Hellolibrary.
  INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES} ${NRRDIO_INCLUDE})

  # Make sure the linker can find the Hello library once it is built.
  LINK_DIRECTORIES(${VISPACK_BINARY_DIR}/lib ${NRRDIO_LIB})

  ADD_EXECUTABLE(linear_test linear_test.cxx)
  TARGET_LINK_LIBRARIES(linear_test matrix util lapack2 lapack1 blas F77 m)

ENDIF(BUILD_BIOMESH3D)
